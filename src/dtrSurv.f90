MODULE INNERS
  IMPLICIT NONE

  PUBLIC

  INTEGER, PARAMETER :: dp = selected_real_kind(15,307)

  INTEGER, SAVE :: ERT ! 0/1 1 = use extremely randomized tree
  INTEGER, SAVE :: minEvent ! minimum number of events in a node
  INTEGER, SAVE :: mTry ! maximum number of covariates to try
  INTEGER, SAVE :: n  ! number of sampled cases
  INTEGER, SAVE :: nAll ! number of cases 
  INTEGER, SAVE :: nLevs ! maximum number of levels 
  INTEGER, SAVE :: nodeSize ! minimum number of cases in a node
  INTEGER, SAVE :: np ! number of covariates
  INTEGER, SAVE :: nrNodes ! maximum number of nodes in a tree
  INTEGER, SAVE :: nt ! number of time points
  INTEGER, SAVE :: nTree ! number of trees
  INTEGER, SAVE :: replace
  INTEGER, SAVE :: rule ! 1 if logrank 0 if truncated mean
  INTEGER, SAVE :: sampleSize
  ! for survival probability, the index of nearest time point
  INTEGER, SAVE :: sIndex
  INTEGER, SAVE :: uniformSplit ! 0/1 1 = random cutoff comes from values

  ! censoring indicator 1 = not censored
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: deltaAll
  ! censoring indicator 1 = not censored
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: delta
 ! number of categories in each np covariate
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nCat

  ! the probability for a random split
  REAL(dp), SAVE :: rs
  ! the fraction above sIndex time that survival time lies
  REAL(dp), SAVE :: sFraction
  ! the stratified random split coefficient
  REAL(dp), SAVE :: stratifiedSplit

  ! time differences
  REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE :: dt
  ! probability mass vector of survival function for sampled cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: pr
  ! probability mass vector of survival function for all cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: prAll 
  ! covariates to be considered for split for sampled cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: x
  ! covariates to be considered for split for all cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: xAll

  ! TRUE = time differences vector has been allocated
  LOGICAL, SAVE :: dtAllocated = .FALSE.
  ! TRUE = all other allocatables have been allocated
  LOGICAL, SAVE :: isAllocated = .FALSE.
  ! TRUE = using survival probability
  LOGICAL, SAVE :: isSurvival = .FALSE.

  TYPE NODE
    INTEGER :: nNode
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: survFunc
    REAL(dp), DIMENSION(:), ALLOCATABLE :: mean
    REAL(dp), DIMENSION(:), ALLOCATABLE :: survProb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: matrix
  END TYPE

  TYPE(NODE), DIMENSION(:), ALLOCATABLE, SAVE :: trees

  TYPE ForestValues
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: survFunc
    REAL(dp), DIMENSION(:), ALLOCATABLE :: mean
    REAL(dp), DIMENSION(:), ALLOCATABLE :: survProb
  END TYPE

  TYPE(ForestValues), SAVE :: forest

  CONTAINS

! sample an array of indices allowing for duplicates
!   nCases: integer, the total number of indices to sample
!   n: integer, the number of indices to draw
! returns an array (1:n) of the indices sampled
FUNCTION sampleWithReplace(nCases, n) RESULT(array)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, INTENT(IN) :: n

  INTEGER, DIMENSION(1:n) :: array

  INTEGER :: i

  REAL(dp) :: rnd

  EXTERNAL :: rnd

  DO i = 1, n
    array(i) = 1 + floor(rnd(0.d0, 1.d0)*nCases)
  END DO

END FUNCTION sampleWithReplace

! sample an array of indices without allowing duplicates
!   nCases: integer, the total number of indices to sample
!   n: integer, the number of indices to draw
! returns an array (1:n) of the indices sampled
FUNCTION sampleWithOutReplace(nCases, n) RESULT(array)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, INTENT(IN) :: n

  INTEGER, DIMENSION(1:n) :: array

  INTEGER :: i, j, cnt

  REAL(dp) :: u

  LOGICAL, DIMENSION(1:nCases) :: used

  REAL(dp) :: rnd

  EXTERNAL :: rnd

  used(:) = .FALSE.

  cnt = 1
  DO WHILE (cnt <= n)
    ! draw a number
    j = 1 + floor(rnd(0.d0, 1.d0)*nCases)
    ! if the number has already been drawn cycle
    if (used(j)) CYCLE
    ! if the number has not previously been drawn, store and set used
    array(cnt) = j
    used(j) = .TRUE.
    cnt = cnt + 1
  END DO

END FUNCTION sampleWithOutReplace


! Identify the optimal split
!   nCases : integer, the number of elements in input casesIn
!   casesIn : integer(:), the indices of the cases in this node
!   nv : integer, the number of covariates to include in search
!   varsIn : integer(:), covariates to include in search
!   splitVar : integer, the index of the selected variable for splitting
!   cutoffBest : real(:), the cutoff (<= go to 'left' node)
!   splitFound : integer, 0 = no split; 1 = found a split
!   casesOut : integer(:), elements of casesIn that go left; ind if yes, 0 
!     otherwise
!   nCuts : integer, the number of cutoff values returned
!   lft : integer, the number of cases in the left node
SUBROUTINE tfindSplit(nCases, casesIn, nv, varsIn, &
                    & splitVar, cutoffBest, splitFound, casesOut, nCuts, lft)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  INTEGER, INTENT(IN) :: nv
  INTEGER, DIMENSION(1:nv), INTENT(IN) :: varsIn
  INTEGER, INTENT(OUT) :: splitVar
  REAL(dp), DIMENSION(1:nLevs), INTENT(OUT) :: cutoffBest
  INTEGER, INTENT(OUT) :: splitFound
  INTEGER, DIMENSION(1:nCases), INTENT(OUT) :: casesOut
  INTEGER, INTENT(OUT) :: nCuts
  INTEGER, INTENT(OUT) :: lft

  INTEGER :: cnt, i, ikv, j, jj, k, kv, l, nUncensored, ptr, rightNode
  INTEGER :: rUnifSet, set, splitLeft, splitLeftFinal, tieCovariate
  INTEGER :: tieValue, variablesTried
  INTEGER, DIMENSION(1:nCases) :: cases, cOut, dSorted, indices, ix, tcases
  INTEGER, DIMENSION(1:nv) :: variables
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind, indSingles, leftCases, rightCases
  INTEGER, DIMENSION(:), ALLOCATABLE :: uncensoredIndices

  REAL(dp) :: cutoff, maxValueSplit, maxValueXm, rUnif, u, valuej
  REAL(dp), DIMENSION(1:nt) :: atRiskj, atRiskLeft, atRiskRight, D, denJ, eventj
  REAL(dp), DIMENSION(1:nt) :: eventsLeft, eventsRight, numJ, pd1, pd2, R, Rcum
  REAL(dp), DIMENSION(1:nCases) :: xSorted
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: prl, prr

  LOGICAL :: randomSplit
  LOGICAL, DIMENSION(:), ALLOCATABLE :: singles

  REAL(dp) :: rnd

  EXTERNAL :: rnd

  ! determine if this is to be a random split
  randomSplit = rnd(0.d0, 1.d0) <= rs

  ! splitFound is flag indicating if a split was found
  ! 0 = no split found; 1 = split found
  splitFound = 0

  !! initialize tracking variables

  ! index of variable with largest critical value
  ! negative number indicates no split found
  splitVar = -1

  ! largest critical value
  maxValueSplit = 0.0

  ! cutoff that gives largest critical value or is randomly selected
  cutoffBest = 0.0

  ! set casesOut to 0 
  casesOut = 0

  ! set number of cutoffs to 0
  nCuts = 0

  ! indices of variables to be explored
  variables = varsIn

  ! location of last available parameter in variables
  ptr = nv

  ! tracks the number of variables tried
  variablesTried = 0

  tcases = (/(i,i=1,nCases)/)

  IF (nCases < 2*nodeSize) RETURN

  DO i = 1, nv

    ! if mTry successful splits already explored, exit
    IF (variablesTried .EQ. mTry) EXIT

    ! randomly select a covariate on which to split
    ikv = 1 + floor(rnd(0.d0, 1.d0)*ptr)
    kv = variables(ikv)

    ! move this index to the end of the list; shift pointer down
    j = variables(ptr)
    variables(ptr) = kv
    variables(ikv) = j

    ptr = ptr - 1

    ! pull appropriate covariate. If un-ordered factor use mean survival time
    ! for each factor level
    IF (nCat(kv) .GT. 1) THEN
      CALL getCovariate(nCases, casesIn, kv, xSorted)
    ELSE
      xSorted = x(casesIn,kv)
    END IF

    cases = casesIn

    ! sort the covariate and track the indices
    CALL qsort4 (xSorted, cases, 1, nCases)
!    CALL hpsort_eps_epw(nCases, xSorted, cases, 1d-8)

    ! sort event indicator data accordingly
    dSorted = delta(cases)

    ! ******************** splitBoundaries ********************
    ! identify minimum cases for left and right splits based on
    ! minimum uncensored cases, minimum node size, and assurance
    ! that all equal valued cases are included in the minimum nodes
    rUnif = 0.d0
    rUnifSet = -1

    ! cases that are not-censored
    uncensoredIndices = pack(tcases, dSorted .EQ. 1)
    nUncensored = size(uncensoredIndices)

    ! if too few cases to meet minimum number of uncensored cases, CYCLE
    IF (nUncensored .LT. (minEvent * 2)) CYCLE

    !! able to split and satisfy minimum number of events in each node

    ! cases to left include all indices up to and including minEvent case
    ! must have at least nodeSize cases
    splitLeft = max(uncensoredIndices(minEvent), nodeSize)

    ! move splitLeft up to include cases with equivalent values of x
    splitLeft = count(xSorted .LE. (xSorted(splitLeft) + 1e-8))

    ! cases to right
    ! include all indices down to and including nUncensored - minEvent + 1 case
    ! must have at least nodeSize cases
    rightNode = min(uncensoredIndices(nUncensored - minEvent + 1), &
                  & nCases - nodeSize + 1)

    ! move rightNode down to include cases with equivalent values of x
    ! splitLeftFinal is the last possible case for the left node
    splitLeftFinal = count(xSorted .LT. xSorted(rightNode))

    ! if the splitLeft index is above the splitLeftFinal index cycle,
    ! split is not possible
    IF (splitLeft .GT. splitLeftFinal) CYCLE

    rUnifSet = 0
    IF ((.NOT. randomSplit) .AND. (ERT .EQ. 1)) THEN

      !************* getUniformSplit *****************

      !! split based on extremely randomized trees and uniform split inputs

      ! if ERT splitLeft = splitLeftFinal, which is the splitting point
      ! and cutoff is set to random value or mid-point depending on uniERT

      !! extremely randomized tree methods used

      IF (uniformSplit .EQ. 0) THEN

        ! if the cutoff is not determined from a uniform distribution
        ! randomly sample available indices to identify the last case
        ! of the left split to define the minimum; cutoff is the
        ! mid-point {x[r] + x[r+1]} / 2

        ! only indices for which x[r] < x[r+1] can be sampled
        ind = (/ (i, i = splitLeft, splitLeftFinal) /)
        singles = xSorted(ind) .LT. xSorted(ind + 1)
        indSingles = pack(ind, singles)
        cnt = size(indSingles)

        splitLeftFinal = indSingles(1 + floor(rnd(0.d0, 1.d0)*cnt))
        ! the last required case in the left node is now splitLeftFinal
        splitLeft = splitLeftFinal
        rUnif = (xSorted(splitLeftFinal) + xSorted(splitLeftFinal+1))/2.0
        rUnifSet = 1

      ELSE IF (uniformSplit .EQ. 1) THEN
        ! randomly select a value in the range of values that satisfy the
        ! allowed cases in the left/right nodes
        rUnif = rnd(0.d0, 1.d0) * (xSorted(splitLeftFinal+1) - &
              & xSorted(splitLeft)) + xSorted(splitLeft)
        rUnifSet = 1

        ! identify the first case that splits to the right of this value
        ! the preceding case is the last case to the left node
        splitLeftFinal = nCases - count(xSorted > rUnif)

        ! the last required case in the left node is now splitLeftFinal
        splitLeft = splitLeftFinal

      END IF

    END IF

    ! -1 is returned if cannot satisfy minimum requirements for nodes
    ! cycle to next covariate
    IF (rUnifSet .EQ. -1) CYCLE

    ! increment the number of covariates that have been explored
    variablesTried = variablesTried + 1

    !***************** maxValue ***************

    ! set initial values for outputs
    set = 0
    maxValueXm = 0.d0
    cutOff = 0.d0

    leftCases = cases(1:(splitLeft-1))
    rightCases = cases(splitLeft:nCases)

    prl = pr(leftCases,:)
    prr = pr(rightCases,:)

    eventsLeft = sum(prl * &
                   & spread(dSorted(1:(splitLeft-1)), 2, nt), DIM = 1)
    eventsRight = sum(prr * &
                    & spread(dSorted(splitLeft:nCases), 2, nt), DIM = 1)

    pd1 = sum(prl, DIM = 1)
    pd2 = sum(prr, DIM = 1)

    atRiskLeft(1) = splitLeft - 1
    atRiskRight(1) = nCases - splitLeft + 1

    DO j = 2, nt
      atRiskLeft(j) = atRiskLeft(j-1) - pd1(j-1)
      atRiskRight(j) = atRiskRight(j-1) - pd2(j-1)
    END DO

    ! if logrank, do calculations that do not depend on node occupancy
    IF (rule == 1) THEN
      CALL logrankSetUp(atRiskLeft, atRiskRight, eventsLeft, eventsRight, &
                      & numJ, denJ)
    END IF


    cnt = 1
    DO j = splitLeft, splitLeftFinal

      ! at risk indicators for jth case
      ! number of events for jth case
      pd1 = prr(cnt,:)
      cnt = cnt + 1
      D = pd1*delta(cases(j))
      Rcum(1) = 0.0
      Rcum(2) = pd1(1)
      DO k = 3, nt
        IF (pd1(k-1) .GT. 1d-8) THEN
          Rcum(k) = Rcum(k-1) + pd1(k-1)
        ELSE 
          Rcum(k) = Rcum(k-1)
        END IF
      END DO

      Rcum = 1.d0 - Rcum 
      ! number at risk

      ! add the jth case to the left node
      atRiskLeft = atRiskLeft + Rcum

      ! remove the jth case from the right node
      atRiskRight = atRiskRight - Rcum

      ! number of events

      ! add the jth case to the left node
      eventsLeft = eventsLeft + D

      ! remove the jth case from the right node
      eventsRight = eventsRight - D

      ! if the case is not the last case with this covariate value, cycle
      IF (xSorted(j) .GE. (xSorted(j+1) - 1d-8)) CYCLE

      ! calculate test statistic
      IF (rule == 1) THEN
        CALL logrank(atRiskLeft, atRiskRight, eventsLeft, eventsRight, numJ, &
                   & denJ, valuej)
      ELSE
        CALL meanSplit(atRiskLeft, atRiskRight, eventsLeft, eventsRight, valuej)
      END IF

      IF ((set .EQ. 0) .OR. (valuej .GT. maxValueXm)) THEN

        ! if first value or value > current max, save
        IF (rUnifSet .EQ. 1) THEN
          cutoff = rUnif
        ELSE
          cutoff = (xSorted(j) + xSorted(j+1))/2.d0
        END IF
        maxValueXm = valuej
        tieValue = 1
        set = 1

      ELSE IF (valuej > (maxValueXm - 1d-8)) THEN
        ! if value is a tie, randomly determine if cutoff should be taken
        tieValue = tieValue + 1

        IF (rnd(0.d0, 1.d0) < (1.d0 / REAL(tieValue))) THEN
          cutoff = (xSorted(j) + xSorted(j+1))/2.d0
        END IF

      END IF

    END DO

    ! if not successful, cycle to next covariate
    ! this condition should never be true
    IF (set .EQ. 0) CYCLE

    ! if successful, determine if it yields the maximum value of the
    ! covariates considered
    IF ((splitVar .EQ. -1) .OR. (maxValueXm .GT. maxValueSplit)) THEN

      ! if first non-zero or largest value, keep cutoff and value and
      ! reset tie counter to 1
      splitVar = kv
      maxValueSplit = maxValueXm
      tieCovariate = 1

      ! count the number of cases in the left node
      lft = count(xSorted .LE. cutoff)
      casesOut = cases

      cutoffBest = 0

      IF (nCat(kv) .LE. 1) THEN
        cutoffBest(1) = cutoff
        nCuts = 1
      ELSE
        ! for factors, identify factor values contained in left cases
        l = 1
        DO j = 1, nCat(kv)
          DO jj = 1, lft
            IF (nint(x(casesOut(jj),kv)) .EQ. j) THEN
              cutoffBest(l) = j
              l = l + 1
              EXIT
            END IF
          END DO
        END DO
        nCuts = l - 1
      END IF

      ! if a random split was triggered, break out of loop over covariates
      IF (randomSplit) EXIT
    
    ELSE IF (maxValueXm .GT. (maxValueSplit-1d-8)) THEN

      ! if equal to current maximum value, increment tie counter and randomly 
      ! select the cutoff with equal probability for each tie
      tieCovariate = tieCovariate + 1

      IF (rnd(0.d0, 1.d0) < (1.0 / tieCovariate)) THEN
        splitVar = kv
        maxValueSplit = maxValueXm
        ! count the number of cases in the left node
        lft = count(xSorted .LE. cutoff)
        casesOut = cases

        cutoffBest = 0

        IF (nCat(kv) .LE. 1) THEN
          cutoffBest(1) = cutoff
          nCuts = 1
        ELSE
          ! for factors, identify factor values contained in left cases
          l = 1
          DO j = 1, nCat(kv)
            DO jj = 1, lft
              IF (nint(x(casesOut(jj),kv)) .EQ. j) THEN
                cutoffBest(l) = j
                l = l + 1
                EXIT
              END IF
            END DO
          END DO
          nCuts = l - 1
        END IF
      END IF
    END IF

  END DO

  ! if no split was possible return
  if (splitVar .EQ. -1) RETURN

  ! if successful at finding a split set flag and return
  splitFound = 1

  RETURN

END SUBROUTINE tfindSplit


! Calculate the Kaplan Meier estimator
! ns integer, the number of time points
! nj real(:), at risk
! oj real(:), events
! z real(:), estimator
SUBROUTINE kaplan(ns, nj, oj, z)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ns
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: nj
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: oj
  REAL(dp), DIMENSION(1:ns), INTENT(OUT) :: z

  INTEGER :: i

  REAL(dp), DIMENSION(1:ns) :: num, on

  num = nj - oj

  z(1) = num(1) / nj(1)

  IF (ns .LT. 2) RETURN

  DO i = 2, ns
    IF (nj(i) > 1d-8) THEN
      z(i) = (num(i)/nj(i)) * z(i-1)
    ELSE
      z(i) = z(i-1)
    END IF
  END DO

END SUBROUTINE

! Truncated Mean test
! N1j: real(:), at risk in group 1
! N2j: real(:), at risk in group 2
! O1j: real(:), events in group 1
! O2j: real(:), events in group 2
! Z: real, truncated mean
SUBROUTINE meanSplit(N1j, N2j, O1j, O2j, Z)
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N2j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O2j
  REAL(dp), INTENT(OUT) :: Z

  INTEGER :: i, m

  REAL(dp), DIMENSION(1:nt) :: E1, E2

  CALL kaplan(nt, N1j, O1j, E1)

  CALL kaplan(nt, N2j, O2j, E2)

  Z = sum((E1 - E2) * dt)
  Z = Z*Z

END SUBROUTINE

! Log rank test set up
! N1j: real(:), at risk in group 1
! N2j: real(:), at risk in group 2
! O1j: real(:), events in group 1
! O2j: real(:), events in group 2
! numJ: real(:), numerator
! denJ: real(:), denominator
SUBROUTINE logRankSetUp(N1j, N2j, O1j, O2j, numJ, denJ)

  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N2j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O2j
  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: numJ
  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: denJ

  INTEGER :: i

  REAL(dp) :: den, Nj, num, Oj, O1
  REAL(dp), DIMENSION(1:nt) :: tdenJ, tNj, tnumJ, tOj

  LOGICAL, DIMENSION(1:nt) :: elg

  numJ = 0.d0
  denJ = 0.d0

  DO i = 1, nt
    IF (N1j(i) .LT. 1d-8) CYCLE
    IF (N2j(i) .LT. 1d-8) CYCLE
    ! time points for which both events have individuals at risk

    ! number of individuals at risk for type 1 or 2 events
    Nj = N1j(i) + N2j(i)

    ! number of events of type 1 or type 2
    Oj = O1j(i) + O2j(i)

    numJ(i) = Oj / Nj

    denJ(i) = numJ(i) * (Nj - Oj) / (Nj * Nj)

  END DO

END SUBROUTINE logrankSetUp

! Log rank test
! N1j: real(:), at risk in group 1
! N2j: real(:), at risk in group 2
! O1j: real(:), events in group 1
! O2j: real(:), events in group 2
! numJ: real(:), numerator
! denJ: real(:), denominator
! Z: real, test value
SUBROUTINE logRank(N1j, N2j, O1j, O2j, numJ, denJ, Z)

  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N2j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O2j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: numJ
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: denJ
  REAL(dp), INTENT(OUT) :: Z

  INTEGER :: i

  REAL(dp) :: den, num

  num = 0.d0
  den = 0.d0
  DO i = 1, nt
    IF (N1j(i) .LT. 1d-8) CYCLE
    IF (N2j(i) .LT. 1d-8) CYCLE
    ! time points for which both events have individuals at risk

    num = num + O1j(i) - numJ(i) * N1j(i)
    den = den + denJ(i) * N1j(i) * N2j(i)

  END DO

  IF (den .GT. 1d-8) THEN
    Z = num * num / den
  ELSE
    Z = 0.d0
  END IF

END SUBROUTINE

! estimate the survival function and mean survival time
!   nCases: integer, the number of elements in casesIn
!   casesIn: integer, the indices of the subset for which the value is 
!     calculated
!   survFunc: real(:), the estimated survival function
!   mean: real, the estimated mean survival 
SUBROUTINE calcValueSingle(nCases, casesIn, survFunc, mean)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: survFunc
  REAL(dp), INTENT(OUT) :: mean

  INTEGER :: i

  REAL(dp), DIMENSION(1:nt) :: Nj, Oj, Rb

  survFunc = 0.d0
  mean = 0.d0

  ! number of at risk cases at each time point
  ! {nt}
  Rb = sum(pr(casesIn,:), DIM = 1)

  Nj(1) = nCases
  DO i = 2, nt
    Nj(i) = Nj(i-1) - Rb(i-1)
  END DO

  ! number of events at each time point
  ! {nt}
  DO i = 1, nt
    Oj(i) = sum(pr(casesIn, i)*delta(casesIn))
  END DO

  ! Kaplan-Meier estimate survival function
  ! {nt}

  CALL kaplan(nt, Nj, Oj, survFunc)

  ! mean survival time 
  mean = sum(survFunc * dt)

  RETURN
END SUBROUTINE

! For factor covariates, calculated the mean survival time and
! use as covariate values for splitting
! nCases integer, number of cases under consideration
! casesIn, integer(:), cases under consideration
! kv, integer, covariate under consideration
! array, real(:), covariate vector as mean survival times
SUBROUTINE getCovariate(nCases, casesIn, kv, array)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  INTEGER, INTENT(IN) :: kv
  REAL(dp), DIMENSION(1:nCases), INTENT(OUT) :: array

  INTEGER :: i
  INTEGER, DIMENSION(1:nCases) :: ix
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind

  REAL(dp) :: mean
  REAL(dp), DIMENSION(1:nt) :: survFunc

  LOGICAL, DIMENSION(1:nCases) :: indices, inSubset

  array = 0.d0

  ! convert covariate to integers to ensure correct equality tests
  ix = nint(x(casesIn,kv))

  DO i = 1, nCat(kv)

    ! identify cases that have this level
    inSubset = ix .EQ. i

    ! ensure that there are individuals in the lth level
    IF (count(inSubset) .EQ. 0) CYCLE

    ! pack indices and estimate survival
    ind = pack(casesIn, inSubset)

    ! calculate the mean survival time for each individual in this subset
    CALL calcValueSingle(size(ind), ind, survFunc, mean)

    WHERE (inSubset) array = mean

  END DO

END SUBROUTINE

! Calculate survival function, mean survival time, and if appropriate 
! survival probability
! nCases, integer, number of cases under consideration
! casesIn, integer(:), indices of cases under consideration
! survFunc, real(:), estimated survival function
! mean, real, estimated mean survival
! survProb, real, estimated survival probability or zero
SUBROUTINE tcalculateValue(nCases, casesIn, survFunc, mean, survProb)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: survFunc
  REAL(dp), INTENT(OUT) :: mean
  REAL(dp), INTENT(OUT) :: survProb

  survFunc = 0.d0
  mean = 0.d0
  survProb = 0.d0

  CALL calcValueSingle(nCases, casesIn, survFunc, mean)

  IF (.NOT. isSurvival) RETURN

  ! estimate survival probability at SurvivalTime
  survProb = survFunc(sIndex) * (1.d0 - sFraction) + &
           & survFunc(sIndex+1) * sFraction
  IF (survProb .LT. 1d-8) survProb = 0.d0

  RETURN
END SUBROUTINE

! grow each tree
! forestSurvFunc, real(:), survival function averaged over forest
! forestMean, real, mean survival averaged over forest
! forestSurvProb, real, survival probability averaged over forest
SUBROUTINE tsurvTree(forestSurvFunc, forestMean, forestSurvProb)
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nAll*nt), INTENT(OUT) :: forestSurvFunc
  REAL(dp), DIMENSION(1:nAll), INTENT(OUT) :: forestMean
  REAL(dp), DIMENSION(1:nAll), INTENT(OUT) :: forestSurvProb

  INTEGER :: i, iTree, j, k, lft, m, nc, ncur, splitFound, splitVar, stat
  INTEGER, DIMENSION(1:sampleSize) :: indices, jdex, xrand
  INTEGER, DIMENSION(1:np) :: newstat, pindices
  INTEGER, DIMENSION(1:np, 1:nrNodes) :: cstat
  INTEGER, DIMENSION(1:nrNodes, 1:2) :: stm
  INTEGER, DIMENSION(1:nAll) :: allStatus
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind, indOut, leftCases, rightCases, pind

  REAL(dp) :: srs
  REAL(dp), DIMENSION(1:nLevs) :: cutoffBest
  REAL(dp), DIMENSION(1:nrNodes) :: mean, survProb
  REAL(dp), DIMENSION(1:nAll) :: xm
  REAL(dp), DIMENSION(1:nt, 1:nrNodes) :: survFunc
  REAL(dp), DIMENSION(1:nrNodes, 1:(5+nLevs)) :: nMatrix
  REAL(dp), DIMENSION(1:nt, 1:nAll) :: tforestSurvFunc

  LOGICAL, DIMENSION(1:np) :: cand
  LOGICAL, DIMENSION(1:nAll) :: sti, tst

  tforestSurvFunc = 0.d0
  forestMean = 0.d0
  forestSurvProb = 0.d0

  DO iTree = 1, nTree

    survFunc = 0.d0
    mean = 0.d0
    survProb = 0.d0
    nMatrix = 0.0
    allStatus = 1

    ! sample data and set local variables x, pr, and delta to the selected
    ! subset
    IF (replace .EQ. 1) THEN
      xrand = sampleWithReplace(nAll, sampleSize)
      n = sampleSize
      x = xAll(xrand,:)
      pr = prAll(xrand,:)
      delta = deltaAll(xrand)
    ELSE IF (nAll .NE. sampleSize) THEN
      xrand = sampleWithoutReplace(nAll, sampleSize)
      n = sampleSize
      x = xAll(xrand,:)
      pr = prAll(xrand,:)
      delta = deltaAll(xrand)
    ELSE
      n = sampleSize
      x = xAll
      pr = prAll
      delta = deltaAll
    END IF

    ! cutoff for identifying covariates to be explored
    srs = stratifiedSplit / REAL(np)

    ! indices for all cases
    indices = (/(i,i=1,n)/)
    jdex = indices

    ! indices for all covariates
    pindices = (/(i,i=1,np)/)

    !! initialize first node

    ! calculate survival function and mean survival of the node
    CALL calcValueSingle(n, indices, survFunc(:,1), mean(1))

    IF (isSurvival) THEN
      ! estimate survival probability at SurvivalTime
      survProb(1) = survFunc(sIndex,1) * (1.d0 - sFraction) + &
                  & survFunc(sIndex+1,1) * sFraction
      IF (survProb(1) .LT. 1d-8) survProb(1) = 0.d0
    END IF 

    ! determine if the node can split based on basic minimum requirements
    if (n .LE. nodeSize .OR. sum(delta(indices)) .LE. 1) THEN
      nMatrix(1,1) = -1
    ELSE
      nMatrix(1,1) = -2
    END IF

    cstat(:,1) = 0

    ! start and finish locations of indices in node
    stm(1,1) = 1
    stm(1,2) = n

    ! location of most recent storage location in matrices/vectors
    ! ncur is incremented when a node successfully splits indicating the
    ! location in the nodes list where base information for the each daughter
    ! is stored
    ncur = 1

    DO k = 1, nrNodes
      ! if k is beyond current node count or 
      ! current node count at limit, break from loop
      IF (k .GT. ncur .OR. ncur .GT. (nrNodes - 2)) EXIT

      ! if node is not to be split, cycle to next node
      IF (nint(nMatrix(k,1)) .EQ. -1) CYCLE

      ! indices for cases contained in node
      ind = jdex(stm(k,1):stm(k,2))

      ! if there are deficient variables, use only these variables
      cand = cstat(:,k) .LT. floor(srs * sum(cStat(:,k)))
      pind = pack(pindices, cand)
      IF (size(pind) .EQ. 0) pind = pindices

      ! split cases
      indOut = ind

      CALL tfindSplit(size(ind), ind, size(pind), pind, splitVar, cutoffBest, &
                    & splitFound, indOut, nc, lft)

      IF (splitFound .EQ. 0 ) THEN
        ! if no split available, set node k as terminal node
        nMatrix(k,1) = -1
        CYCLE
      END IF

      ! set node k to be interior (i.e. has split) 
      nMatrix(k,1) = -3

      ! add split information to node
      nMatrix(k,4) = pindices(splitVar)
      nMatrix(k,5) = nc
      nMatrix(k,6:(6+nc-1)) = cutoffBest(1:nc)
      nMatrix(k,2) = ncur + 1
      nMatrix(k,3) = ncur + 2

      ! increment the times the variable was used in a split
      newstat = cStat(:,k)
      newstat(nint(nMatrix(k,4))) = newstat(nint(nMatrix(k,4))) + 1

      ! store new case order in jdex
      jdex(stm(k,1):stm(k,2)) = indOut

      !! left node 

      ncur = ncur + 1

      ! index boundaries for cases in left node
      stm(ncur,1) = stm(k,1)
      stm(ncur,2) = stm(k,1) + lft - 1

      leftCases = jdex(stm(ncur,1):stm(ncur,2))

      ! get basic node information for left daughter
      CALL calcValueSingle(size(leftCases), leftCases, survFunc(:,ncur), &
                         & mean(ncur))

      IF (isSurvival) THEN
        ! estimate survival probability at SurvivalTime
        survProb(ncur) = survFunc(sIndex,ncur) * (1.d0 - sFraction) + &
                       & survFunc(sIndex+1,ncur) * sFraction
        IF (survProb(ncur) .LT. 1d-8) survProb(1) = 0.d0
      END IF 

      IF (size(leftCases) .LE. nodeSize .OR. sum(delta(leftCases)) .LE. 1) THEN
        ! if the number of cases in the node is at or below the minimum required
        ! or the number of uncensored event is only 1
        ! status is terminal
        nMatrix(ncur,1) = -1
      ELSE
        nMatrix(ncur,1) = -2
      END IF

      cstat(:,ncur) = newstat

      !! right node 

      ncur = ncur + 1

      ! index boundaries for cases in right node
      stm(ncur,1) = stm(k,1) + lft
      stm(ncur,2) = stm(k,2)

      ! retrieve left and right cases
      rightCases = jdex(stm(ncur,1):stm(ncur,2))

      ! calculate survival function and mean survival time for right node
      CALL calcValueSingle(size(rightCases), rightCases, survFunc(:,ncur), &
                         & mean(ncur))

      IF (isSurvival) THEN
        survProb(ncur) = survFunc(sIndex,ncur) * (1.d0 - sFraction) + &
                       & survFunc(sIndex+1,ncur) * sFraction
        IF (survProb(ncur) .LT. 1d-8) survProb(1) = 0.d0
      END IF 

      IF (size(rightCases) .LE. nodeSize .OR. &
        & sum(delta(rightCases)) .LE. 1) THEN
        ! if the number of cases in the node is at or below the minimum required
        ! or the number of uncensored event is only 1
        ! status is terminal
        nMatrix(ncur,1) = -1
      ELSE
        nMatrix(ncur,1) = -2
      END IF

      cstat(:,ncur) = newstat

      ! retrieve the variable on which data is split
      m = nint(nMatrix(k,4))

      ! retrieve the covariate
      xm = xAll(:,m)

      IF (nCat(m) .LE. 1) THEN
        ! if a numeric variable, use the cutoff value to
        ! identify if individual i goes left
        tst = xm .LE. nMatrix(k,6)
      ELSE
        ! if an unordered factor, use category to identify if individual
        ! i goes left
        tst = .FALSE.
        DO j = 1, nint(nMatrix(k,5))
          tst = tst .OR. nint(xm) .EQ. nint(nMatrix(k,j+5))
        END DO
      END IF

      DO j = 1, nAll
        IF (allStatus(j) .NE. k) CYCLE
        IF (tst(j)) THEN
          allStatus(j) = nint(nMatrix(k,2))
        ELSE 
          allStatus(j) = nint(nMatrix(k,3))
        END IF
      END DO

    END DO

    ! ensure that all nodes that are not "interior" are "terminal"
    WHERE (nMatrix(:,1) .EQ. -2) nMatrix(:,1) = -1

    trees(iTree)%survFunc = survFunc(:,1:ncur)
    trees(iTree)%mean = mean(1:ncur)
    trees(iTree)%survProb = survProb(1:ncur)
    trees(iTree)%matrix = nMatrix(1:ncur,:)
    trees(iTree)%nNode = ncur

    tforestSurvFunc = tforestSurvFunc + survFunc(:,allStatus)
    forestMean = forestMean + mean(allStatus)
    forestSurvProb = forestSurvProb + survProb(allStatus)

  END DO

  forestSurvFunc = reshape(tforestSurvFunc, (/nt*nAll/)) / nTree
  forestMean = forestMean / nTree
  forestSurvProb = forestSurvProb / nTree

END SUBROUTINE tSurvTree

! use tree structure to predict for newData
SUBROUTINE predict(iTree, nr, nc, newData)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iTree
  INTEGER, INTENT(IN) :: nr
  INTEGER, INTENT(IN) :: nc
  REAL(dp), DIMENSION(:,:) :: newData

  INTEGER :: i, j, m
  INTEGER, DIMENSION(1:nr) :: ll, stat

  REAL(dp), DIMENSION(1:nr) :: xm
  REAL(dp), DIMENSION(1:nr,1:nc) :: nData

  LOGICAL, DIMENSION(1:nr) :: sti, tst

  nData = reshape(newData, (/nr,nc/))

  ! column 1 is 0/1 indicating interior/terminal
  ! column 2 is the index of left
  ! column 3 is the index of right
  ! column 4 is the split variable
  ! column 5 is the number of cutoffs
  ! column 6:nCol are the cutoff values

  ! begin with every case in node 1
  stat = 1

  DO i = 1, trees(iTree)%nNode

    ! if terminal cycle
    IF (nint(trees(iTree)%matrix(i,1)) .EQ. -1) CYCLE

    ! identify individuals in this node
    sti = stat .EQ. i

    ! retrieve the variable on which data is split
    m = nint(trees(iTree)%matrix(i,4))

    ! retrieve the covariate
    xm = nData(:,m)

    IF (nCat(m) .LE. 1) THEN
      ! if a numeric variable, use the cutoff value to
      ! identify if individual i goes left
      tst = xm .LE. trees(iTree)%matrix(i,6)
    ELSE
      ! if an unordered factor, use category to identify if individual
      ! i goes left
      tst = .FALSE.
      DO j = 1, nint(trees(iTree)%matrix(i,5))
        tst = tst .OR. nint(xm) .EQ. nint(trees(iTree)%matrix(i,j+5))
      END DO
    END IF

    WHERE (tst .AND. sti) stat = nint(trees(iTree)%matrix(i,2))
    WHERE ( (.NOT. tst) .AND. sti) stat = nint(trees(iTree)%matrix(i,3))

  END DO

  forest%survFunc = forest%survFunc + trees(iTree)%survFunc(:,stat)
  forest%mean = forest%mean + trees(iTree)%mean(stat)
  forest%survProb = forest%survProb + trees(iTree)%survProb(stat)

  RETURN

END SUBROUTINE


END MODULE INNERS

! n, integer, number of cases in data
! np, integer, number of covariates in data
! xt, integer(:), covariates
! nCat, integer(:), number of levels in each covariate (0 = continuous,
!   1 = ordered factors, 2+ = unordered factor)
! nt, integer, number of time points
! nNodes, integer, number of nodes
! tsurvFunc, real(:), survival functions for each node
! mean, real(:), mean survival of each node
! survProb, real(:), survival probability of each node
! nCols, integer, number of columns in node matrix
! tnodes, real(:), node information
! predSurvFunc, real(:), predicted survival functions
! predMean, real(:), predicted mean survival
! predSurvProb, real(:), predicted survival probabilities
SUBROUTINE predictSurvTree(n, np, xt, nCat, nt, nNodes, tsurvFunc, mean, &
  & survProb, nCols, tnodes, predSurvFunc, predMean, predSurvProb)
  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = selected_real_kind(15,307)

  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: np
  REAL(dp), DIMENSION(1:n*np), INTENT(IN) :: xt
  INTEGER, DIMENSION(1:np), INTENT(IN) :: nCat
  INTEGER, INTENT(IN) :: nt
  INTEGER, INTENT(IN) :: nNodes
  REAL(dp), DIMENSION(1:nt*nNodes), INTENT(IN) :: tsurvFunc
  REAL(dp), DIMENSION(1:nNodes), INTENT(IN) :: mean
  REAL(dp), DIMENSION(1:nNodes), INTENT(IN) :: survProb
  INTEGER, INTENT(IN) :: nCols
  REAL(dp), DIMENSION(1:nNodes*nCols), INTENT(IN) :: tnodes
  REAL(dp), DIMENSION(1:n*nt), INTENT(OUT) :: predSurvFunc
  REAL(dp), DIMENSION(1:n), INTENT(OUT) :: predMean
  REAL(dp), DIMENSION(1:n), INTENT(OUT) :: predSurvProb

  INTEGER :: areTerminal, i, j, m
  INTEGER, DIMENSION(1:n) :: ll, stat

  REAL(dp), DIMENSION(1:n) :: xm
  REAL(dp), DIMENSION(1:nt, 1:n) :: tsurv
  REAL(dp), DIMENSION(1:n, 1:np) :: x
  REAL(dp), DIMENSION(1:nNodes, 1:nCols) :: nodes
  REAL(dp), DIMENSION(1:nt, 1:nNodes) :: survFunc

  LOGICAL, DIMENSION(1:n) :: sti, tst

  x = reshape(xt, (/n, np/))
  nodes = reshape(tnodes, (/nNodes, nCols/))
  survFunc = reshape(tsurvFunc, (/nt, nNodes/))

  tsurv = 0.d0

  ! column 1 is indicator of interior/terminal
  ! column 2 is the index of left
  ! column 3 is the index of right
  ! column 4 is the split variable
  ! column 5 is the number of cutoffs
  ! column 6:nCol are the cutoff values

  ! begin with every case in node 1
  ! stat will contain the terminal node to which each case belongs
  stat = 1

  DO i = 1, nNodes
    ! if terminal cycle

    IF (nint(nodes(i,1)) .EQ. -1) CYCLE

    ! identify individuals in this node
    sti = stat .EQ. i

    ! retrieve the variable on which data is split
    m = nint(nodes(i,4))

    ! retrieve the covariate
    xm = x(:,m)

    IF (nCat(m) .LE. 1) THEN
      ! if a numeric variable, use the cutoff value to
      ! identify if individual i goes left
      tst = xm .LE. nodes(i,6)
    ELSE
      ! if an unordered factor, use category to identify if individual
      ! i goes left
      tst = .FALSE.
      DO j = 1, nint(nodes(i,5))
        tst = tst .OR. nint(xm) .EQ. nint(nodes(i,j+5))
      END DO
    END IF

    WHERE (tst .AND. sti) stat = nint(nodes(i,2))
    WHERE ( (.NOT. tst) .AND. sti) stat = nint(nodes(i,3))

  END DO

  ! retrieve appropriate values based on terminal node
  tsurv = survFunc(:,stat)
  predSurvFunc = reshape(tsurv,(/n*nt/))
  predMean = mean(stat)
  predSurvProb = survProb(stat)

END SUBROUTINE

! set up basic information for the module
! t_nt, integer, the number of time points
! t_dt, real(:), the time differences between time points
! t_rs, real, the probability for a random split
! t_ERT, integer, the indicator of extremely randomized trees
! t_uniformSplit, integer, the indicator of method for determining cut-off
!   when using ERT
! t_nodeSize, integer, the minimum number of cases in each node
! t_minEvent, integer, the minimum number of events in each node
! t_rule, integer, 0 = mean, 1 = logrank
! t_sIndex, integer, the indices of time points that is closest to the 
!   requested survival time
! t_sFraction, real, the fractional distance between time points the the 
!   requested survival time
! t_stratifiedSplit, real, the coefficient for determining stratification
! t_replace, integer, indicator of sampling with replacement
SUBROUTINE setUpBasics(t_nt, t_dt, t_rs, t_ERT, t_uniformSplit, t_nodeSize, &
                     & t_minEvent, t_rule, t_sIndex, t_sFraction, &
                     & t_stratifiedSplit, t_replace)

  USE INNERS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: t_nt
  REAL(dp), DIMENSION(1:t_nt), INTENT(IN) :: t_dt
  REAL(dp), INTENT(IN) :: t_rs
  INTEGER, INTENT(IN) :: t_ERT
  INTEGER, INTENT(IN) :: t_uniformSplit
  INTEGER, INTENT(IN) :: t_nodeSize
  INTEGER, INTENT(IN) :: t_minEvent
  INTEGER, INTENT(IN) :: t_rule
  INTEGER, INTENT(IN) :: t_sIndex
  REAL(dp), INTENT(IN) :: t_sFraction
  REAL(dp), INTENT(IN) :: t_stratifiedSplit
  INTEGER, INTENT(IN) :: t_replace

  nt = t_nt
  sIndex = t_sIndex
  sFraction = t_sFraction

  isSurvival = sIndex > 0

  IF (dtAllocated) DEALLOCATE(dt)

  ALLOCATE(dt(1:nt))

  dtAllocated = .TRUE.

  dt = t_dt

  rs = t_rs
  ERT = t_ERT
  uniformSplit = t_uniformSplit
  nodeSize = t_nodeSize
  minEvent = t_minEvent
  rule = t_rule
  stratifiedSplit = t_stratifiedSplit
  replace = t_replace

END SUBROUTINE setUpBasics

! set up basic information for the module that is step dependent
! t_n, integer, the number of cases under consideration
! t_np, integer, the number of covariates
! t_x, real(:), the covariates
! t_pr, real(:), the probability mass vector of survival function
! t_delta, integer(:), the indicator of censoring
! t_mTry, integer, the maximum number of covariates to try for splitting
! t_nCat, integer(:), the number of categories in each covariate
! t_sampleSize, integer, the number of cases to sample for each tree
! t_ntree, integer, the number of trees in the forest
! t_nrNodes, integer, the maximum number of nodes
SUBROUTINE setUpInners(t_n, t_np, t_x, t_pr, t_delta, t_mTry, t_nCat, &
                     & t_sampleSize, t_nTree, t_nrNodes)

  USE INNERS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: t_n
  INTEGER, INTENT(IN) :: t_np
  REAL(dp), DIMENSION(1:t_n*t_np), INTENT(IN) :: t_x
  REAL(dp), DIMENSION(1:nt*t_n), INTENT(IN) :: t_pr
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_delta
  INTEGER, INTENT(IN) :: t_mTry
  INTEGER, DIMENSION(1:t_np), INTENT(IN) :: t_nCat
  INTEGER, INTENT(IN) :: t_sampleSize
  INTEGER, INTENT(IN) :: t_nTree
  INTEGER, INTENT(IN) :: t_nrNodes

  nAll = t_n
  np = t_np

  IF (isAllocated) THEN
    DEALLOCATE(xAll, prAll, deltaAll, nCat, forest%survFunc, forest%mean,  &
             & forest%survProb, trees)
  END IF

  ALLOCATE(xAll(1:nAll,1:np))
  ALLOCATE(prAll(1:nAll, 1:nt))
  ALLOCATE(deltaAll(1:nAll))
  ALLOCATE(nCat(1:np))

  isAllocated = .TRUE.

  xAll = reshape(t_x, (/nAll,np/))
  prAll = reshape(t_pr, (/nAll,nt/))
  deltaAll = t_delta
  nCat = t_nCat
  nLevs = max(maxval(nCat),1)

  mTry = t_mTry
  sampleSize = t_sampleSize

  ALLOCATE(forest%survFunc(1:nt, 1:nAll))
  ALLOCATE(forest%mean(1:nAll))
  ALLOCATE(forest%survProb(1:nAll))
  forest%survFunc = 0.d0
  forest%mean = 0.d0
  forest%survProb = 0.d0

  nTree = t_nTree

  ALLOCATE(trees(1:nTree))

  nrNodes = t_nrNodes

END SUBROUTINE setUpInners

! access function for calculating forest
SUBROUTINE survTree(tSurvFunc, mean, survProb)
  USE INNERS
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nrNodes*nt), INTENT(OUT) :: tSurvFunc
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: mean
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: survProb

  CALL tsurvTree(tSurvFunc, mean, survProb)

END SUBROUTINE

! retrieve dimensions of node matrix for the iTree-th tree
SUBROUTINE treeDim(iTree, nr, nc)
  USE INNERS
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iTree
  INTEGER, INTENT(OUT) :: nr
  INTEGER, INTENT(OUT) :: nc

  nr = size(trees(iTree)%matrix,1)
  nc = size(trees(iTree)%matrix,2)

END SUBROUTINE

! retrieve the nodes, survival function, mean survival, and survival probability
! for the iTree-th tree
SUBROUTINE getTree(iTree, nr, nc, nodes, survFunc, mean, survProb)
  USE INNERS
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iTree
  INTEGER, INTENT(IN) :: nr
  INTEGER, INTENT(IN) :: nc
  REAL(dp), DIMENSION(1:nr*nc), INTENT(OUT) :: nodes
  REAL(dp), DIMENSION(1:nt*nr), INTENT(OUT) :: survFunc
  REAL(dp), DIMENSION(1:nr), INTENT(OUT) :: mean
  REAL(dp), DIMENSION(1:nr), INTENT(OUT) :: survProb

  nodes = reshape(trees(iTree)%matrix, (/nr*nc/))
  survFunc = reshape(trees(iTree)%survFunc, (/nt*nr/))
  mean = trees(iTree)%mean
  survProb = trees(iTree)%survProb

END SUBROUTINE

!FUNCTION RND(a, b) RESULT(dRes)
!  implicit none

!  REAL(dp), INTENT(IN) :: a
!  REAL(dp), INTENT(IN) :: b
!  REAL(dp) :: dRes

!  CALL random_number(dRes)
!END FUNCTION
