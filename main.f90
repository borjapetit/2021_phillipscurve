
! ******************************************************************************
! THIS CODE CONTROLS THE EXECUTION OF THE PROGRAM AND CONTAINS THE SUBROUTINES
! TO SOLVE THE STEADY STATE, THE DYNAMICS AND TO CALIBRATE THE MODEL.
! ******************************************************************************

! Compilation command:
!   cd ..
!   gfortran -fopenmp -O3 -ffixed-line-length-150 -fmax-stack-var-size=1000000 -J $(pwd)/compiledfiles toolkit.f90 parameters.f90 solution.f90 dynamics.f90 main.f90 -o lpw

PROGRAM main
USE omp_lib
USE parameters
IMPLICIT NONE
INTEGER           :: i,j,solmet,ip,is,iw,iz,solmet2,cnoise(4),cinfl(5)
REAL(rp)          :: time1,time0,timemp1,timemp0,kappa_pi0,kappa_w0
REAL(rp)          :: PRICES1(3),ZEROS1(3),PARSMAX(7),PARSMIN(7)
CHARACTER(LEN=8)  :: date
CHARACTER(LEN=10) :: time
CHARACTER(LEN=1)  :: vers0,vers1

! ******************************************************************************
! INITIALIZE COMPUTATION

! Time and date of execution
CALL DATE_AND_TIME(DATE=date,TIME=time)

! Set working directory
CALL CHDIR(TRIM(ADJUSTL(path))//"textfiles/")

! Start timing
CALL CPU_TIME(time0)

! Start timing (for OpenMp calculations)
timemp0    = omp_get_wtime()

! Set number of threads for parallel computation
numthreads = omp_get_max_threads( )

! ******************************************************************************
! PRINT HEADER

PRINT *, '  '
PRINT *, ' **************************************************************** '
PRINT *, ' MONETARY POLICY IMPLICATIONS OF STATE-DEPENDENT PRICES AND WAGES '
PRINT *, ' Anton Nakov, James Costain and Borja Petit                       '
PRINT *, ' 2020                                                             '
PRINT *, ' **************************************************************** '
PRINT *, '  '
PRINT *, ' Date: ',date(7:8),'/',date(5:6),'/',date(1:4)
PRINT *, ' Time: ',time(1:2),':',time(3:4)
PRINT *, '  '
WRITE(*,'(A,I3)') '  Threads: ',numthreads
PRINT *, '  '
PRINT *, ' **************************************************************** '

! ******************************************************************************
! INITIALIZE PARAMETERS

! Read calibrated parameters
OPEN(unit=1,file="calibparams.txt",action='read')
  READ(1,*)  lbar
  READ(1,*)  rhobar
  READ(1,*)  kappa_pi
  READ(1,*)  kappa_lambda
  READ(1,*)  kappa_w
  READ(1,*)  kappa_rho
  READ(1,*)  rho_z
  READ(1,*)  stdMC_z
  READ(1,*)  rho_s
  READ(1,*)  stdMC_s
CLOSE (1)

! Baseline adjustment cost parameters (for experiments)
kappa_pi0 = kappa_pi
kappa_w0  = kappa_w

! Extra inflation rate (for experiements)
mu0 = zero

! Build grids and transition matrices
CALL SET_MATS( )

! Fill vector with empirical moments from the data
CALL SET_MOMENTS( )

! ******************************************************************************
! DEFINE WHAT TO DO AND INITIALIZE SOLUTION

vers0 = " "   ! Model version - adjustment costs
vers1 = " "   ! Model version - inflation rate

! Ask what to do
WRITE(*,'(A)') '                                                 '
WRITE(*,'(A)') '  What do you want to do?                        '
WRITE(*,'(A)') '                                                 '
WRITE(*,'(A)') '    (1) Solve steady state                       '
WRITE(*,'(A)') '    (2) Solve steady state and dynamics          '
WRITE(*,'(A)') '    (3) Calibrate                                '
WRITE(*,'(A)') '    (4) Solve for diff. noise parameters         '
WRITE(*,'(A)') '    (5) Solve for diff. inflation rates          '
WRITE(*,'(A)') '    (6) Solve for all cases                      '
WRITE(*,'(A)') '    (7) Solve for pre and post inflation rates   '
WRITE(*,'(A)') '                                                 '
WRITE(*,'(A)',ADVANCE="NO") '    --> Your choice (1-7): ' ; READ (*,*) solmet

! Ask again if incorrect choice
IF (solmet.LT.1 .AND. solmet.GT.7) THEN
  WRITE(*,'(A)',ADVANCE="NO") '    Incorrect choice: from 1 to 7. ' ; GOTO 9
END IF

! If solving steady state and/or dynamics, choose adjustment costs and inflation
IF (solmet.EQ.1 .OR. solmet.EQ.2) THEN

  ! Choose adjustment costs
  vers0 = " "
  WRITE(*,'(A)',ADVANCE="YES") '                                              '
  WRITE(*,'(A)',ADVANCE="YES") '  Which version?                              '
  WRITE(*,'(A)',ADVANCE="YES") '                                              '
  WRITE(*,'(A)',ADVANCE="YES") '    (1) Sticky prices and sticky wages        '
  WRITE(*,'(A)',ADVANCE="YES") '    (2) Semi-flexible prices and sticky wages '
  WRITE(*,'(A)',ADVANCE="YES") '    (3) Flexible prices and sticky wages      '
  WRITE(*,'(A)',ADVANCE="YES") '    (4) Sticky prices and sticky wages        '
  WRITE(*,'(A)',ADVANCE="YES") '    (5) Sticky prices and semi-lexible wages  '
  WRITE(*,'(A)',ADVANCE="YES") '    (6) Flexible prices and flexible wages    '
  WRITE(*,'(A)',ADVANCE="YES") '                                              '
  WRITE(*,'(A)',ADVANCE="NO")  '    --> Your choice (1-6): ' ; READ (*,*) j
  IF (j.LT.1 .AND. j.GT.6) THEN
    WRITE (*,'(A)',ADVANCE="NO") '    Incorrect version: from 1 to 6. '
    GOTO 9
  END IF
  WRITE(vers0,'(I1)') j
  CALL SETKAPPAS(vers0,kappa_pi0,kappa_w0)

  ! Choose inflation rate
  vers1 = " "
  WRITE(*,'(A)',ADVANCE="YES") '                                    '
  WRITE(*,'(A)',ADVANCE="YES") '  Which inflation rate?             '
  WRITE(*,'(A)',ADVANCE="YES") '                                    '
  WRITE(*,'(A)',ADVANCE="YES") '    (0) Inflation = 2% (baseline)   '
  WRITE(*,'(A)',ADVANCE="YES") '    (1) Inflation = -1%             '
  WRITE(*,'(A)',ADVANCE="YES") '    (2) Inflation = 0%              '
  WRITE(*,'(A)',ADVANCE="YES") '    (3) Inflation = 4%              '
  WRITE(*,'(A)',ADVANCE="YES") '    (4) Inflation = 8%              '
  WRITE(*,'(A)',ADVANCE="YES") '    (5) Inflation = -2%             '
  WRITE(*,'(A)',ADVANCE="YES") '    (6) Inflation = 1%              '
  WRITE(*,'(A)',ADVANCE="YES") '                                    '
  WRITE(*,'(A)',ADVANCE="NO" ) '    --> Your choice (0-6): ' ; READ (*,*) j
  IF (j.LT.0 .OR. j.GT.6) THEN
    WRITE (*,'(A)',ADVANCE="NO") '    Incorrect version: from 0 to 6. ' ; GOTO 9
  END IF
  WRITE(vers1,'(I1)') j
  CALL SETINFLATION(vers1)

END IF

IF (solmet.EQ.4 .OR. solmet.EQ.5 .OR. solmet.EQ.6) THEN

  WRITE(*,'(A)',ADVANCE="YES") '  What to solve?                   '
  WRITE(*,'(A)',ADVANCE="YES") '                                   '
  WRITE(*,'(A)',ADVANCE="YES") '    (1) Steady-state               '
  WRITE(*,'(A)',ADVANCE="YES") '    (2) Steady-state and dynamics  '
  WRITE(*,'(A)',ADVANCE="YES") '                                   '
  WRITE(*,'(A)',ADVANCE="NO")  '    --> Your choice (1-7): ' ; READ (*,*) solmet2
  IF (j.LT.1 .AND. j.GT.2) THEN
    WRITE (*,'(A)',ADVANCE="NO") '    Incorrect choice: either 1 or 2. '
    GOTO 9
  END IF

END IF


WRITE(*,'(A)') '  '

PRINT * , ' **************************************************************** '
PRINT * , '  '

! ******************************************************************************
! IMPLEMENT DESIRED SOLUTION

! Solve steady state for given version
IF (solmet.eq.1) THEN

  CALL PRINTVERSION("V"//vers0//vers1)
  CALL SOLVESTEADY("V"//vers0//vers1,2)

! Compute Jacobian for given version
ELSE IF (solmet.eq.2) THEN

  CALL PRINTVERSION("V"//vers0//vers1)
  CALL SOLVEDYN("V"//vers0//vers1)

! Calibrate the parameters
ELSE IF (solmet.eq.3) THEN

  CALL PRINTVERSION("VC0")
  CALL CALIBRATE( )

! Compute the steady state and the jacobian for the different values of the noise parameters
ELSEIF (solmet.EQ.4) THEN ; WRITE(vers1,'(I1)') 0
  DO j=1,6 ; WRITE(vers0,'(I1)') j
    CALL SETKAPPAS(vers0,kappa_pi0,kappa_w0)
    CALL SETINFLATION(vers1)
    CALL PRINTVERSION("V"//vers0//vers1)
    IF (solmet2.eq.1) CALL SOLVESTEADY("V"//vers0//vers1)
    IF (solmet2.eq.2) CALL SOLVEDYN("V"//vers0//vers1)
  END DO
  GOTO 9

! Compute the steady state and the jacobian for different values of the inflation rate
ELSEIF (solmet.EQ.5) THEN ; WRITE(vers0,'(I1)') 1
  DO j = 0,6 ; WRITE(vers1,'(I1)') j
    CALL SETKAPPAS(vers0,kappa_pi0,kappa_w0)
    CALL SETINFLATION(vers1)
    CALL PRINTVERSION("V"//vers0//vers1)
    IF (solmet2.eq.1) CALL SOLVESTEADY("V"//vers0//vers1)
    IF (solmet2.eq.2) CALL SOLVEDYN("V"//vers0//vers1)
  END DO
  GOTO 9

! Compute the steady state for all cases
ELSEIF (solmet.EQ.6) THEN ;
  cnoise = (/ 1, 3, 5, 6    /)
  cinfl  = (/ 5, 2, 0, 3, 4 /)
  DO i = 1,4 ; WRITE(vers0,'(I1)') cnoise(i)
  DO j = 1,5 ; WRITE(vers1,'(I1)') cinfl(j)
  !DO i = 1,6 ; WRITE(vers0,'(I1)') i
  !DO j = 0,6 ; WRITE(vers1,'(I1)') j
    CALL SETKAPPAS(vers0,kappa_pi0,kappa_w0)
    CALL SETINFLATION(vers1)
    CALL PRINTVERSION("V"//vers0//vers1)
    IF (solmet2.eq.1) CALL SOLVESTEADY("V"//vers0//vers1)
    IF (solmet2.eq.2) CALL SOLVEDYN("V"//vers0//vers1)
  END DO
  END DO
  GOTO 9

! Compute the steady state and the jacobian for cases 7 and 8 (pre and post 2000)
ELSEIF (solmet.EQ.7) THEN
  vers0 = "1" ! vers0 = 1: baseline cost parameters
  CALL SETKAPPAS(vers0,kappa_pi0,kappa_w0)
  DO j=7,8 ; WRITE(vers1,'(I1)') j
    CALL SETINFLATION(vers1)
    CALL PRINTVERSION("V"//vers0//vers1)
    CALL SOLVEDYN("V"//vers0//vers1)
  END DO
  GOTO 9

END IF

! ******************************************************************************
! PRINTING STATISTICS

PRINT * , '  '
PRINT *, ' **************************************************************** '
PRINT * , '  '
PRINT ('(A2,A7,F12.10)'), '  ' , 'cbar   ' , cbar
PRINT ('(A2,A7,F12.10)'), '  ' , 'nbar   ' , nbar
PRINT ('(A2,A7,F12.10)'), '  ' , 'wbar   ' , wbar
PRINT * , '  '
PRINT ('(A2,60(A1))'     ) , '  ', ('-',j=1,46)
PRINT ('(A2,A15,A14,A3,A14)' ) , '  ', '  ','      Prices  ', '  |', '      Wages   '
PRINT ('(A2,A15,2(A7,A7,A3))') , '  ', ' ','   Data','  Model', '  |', '   Data',' Model'
PRINT ('(A2,60(A1))'     ) , '  ', ('-',j=1,46)
DO j=1,9
  PRINT ('(A2,A15,2(F7.2,F7.2,A3))') , '  ' , MOMPRINTNAME(j) , &
  MOMPRINTDATA(j)   , MOMPRINTMODEL(j) , '  |' , &
  MOMPRINTDATA(j+19) , MOMPRINTMODEL(j+19)
END DO
PRINT ('(A2,60(A1))'     ) , '  ', ('-',j=1,46)
DO j=10,13
  PRINT ('(A2,A15,2(A7,F7.2,A3))') , '  ' , MOMPRINTNAME(j) , &
  '   -- ' , MOMPRINTMODEL(j) , '  |' , &
  '   -- ' , MOMPRINTMODEL(j+19)
END DO
PRINT ('(A2,60(A1))'     ) , '  ', ('-',j=1,46)
PRINT ('(A2,A15,2(A7,F7.2,A3))') , '  ' , MOMPRINTNAME(14) , &
  '   -- ' , MOMPRINTMODEL(14) , '  |'
PRINT ('(A2,60(A1))'     ) , '  ', ('-',j=1,32)
PRINT * , '  '
PRINT *, ' **************************************************************** '
PRINT * , '  '

! ******************************************************************************
! PRINTING TIMING:

9 CALL CPU_TIME(time1)
  timemp1 = omp_get_wtime()

PRINT * , '            '
PRINT * , ' Completed! '
PRINT * , '            '
IF (time1   - time0   .lt. DBLE(60.0)) PRINT * , ' Time:' , time1-time0                  ,'segs'
IF (time1   - time0   .ge. DBLE(60.0)) PRINT * , ' Time:' , (time1-time0)/DBLE(60.0)     ,'mins'
IF (timemp1 - timemp0 .lt. DBLE(60.0)) PRINT * , ' Time:' , timemp1-timemp0              ,'segs'
IF (timemp1 - timemp0 .ge. DBLE(60.0)) PRINT * , ' Time:' , (timemp1-timemp0)/DBLE(60.0) ,'mins'
PRINT * , '            '
PRINT * , '            '
PRINT * , '            '

CONTAINS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SE PARAMETER VALUES FOR EXPERIMENTS
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! SET DIFFERENT VALUES FOR THE COST OF DECISION-MAKING
SUBROUTINE SETKAPPAS(vers,kp0,kw0)
  IMPLICIT NONE
  REAL(rp)         , INTENT(IN) :: kp0,kw0
  CHARACTER(LEN=1) , INTENT(IN) :: vers
  IF (vers.EQ."1") THEN
    kappa_pi = kp0
    kappa_w  = kw0
  ELSE IF (vers.EQ."2") THEN
    kappa_pi = kp0/diez
    kappa_w  = kw0
  ELSE IF (vers.EQ."3") THEN
    kappa_pi = kp0/cien
    kappa_w  = kw0
  ELSE IF (vers.EQ."4") THEN
    kappa_pi = kp0
    kappa_w  = kw0/diez
  ELSE IF (vers.EQ."5") THEN
    kappa_pi = kp0
    kappa_w  = kw0/cien
  ELSE IF (vers.EQ."6") THEN
    kappa_pi = kp0/cien
    kappa_w  = kw0/cien
  END IF
  kappa_lambda = kappa_pi
  kappa_rho    = kappa_w
  RETURN
END SUBROUTINE SETKAPPAS

! SET DIFFERENT VALUES FOR THE STEADY-STATE INFLATION RATE
SUBROUTINE SETINFLATION(mucase)
  IMPLICIT NONE
  CHARACTER(LEN=1) , INTENT(IN) :: mucase
  IF (mucase.EQ." ") mu0 = zero                                                   ! Inflation rate ~  2%   ---> 0% annual more
  IF (mucase.EQ."0") mu0 = zero                                                   ! Inflation rate ~  2%   ---> 0% annual more
  IF (mucase.EQ."1") mu0 = (((mu**DBLE(12.0))-DBLE(0.03))**(one/DBLE(12.0))) - mu ! Inflation rate ~ -1%   ---> 3% annual less
  IF (mucase.EQ."2") mu0 = (((mu**DBLE(12.0))-DBLE(0.02))**(one/DBLE(12.0))) - mu ! Inflation rate ~ +0%   ---> 2% annual less
  IF (mucase.EQ."3") mu0 = (((mu**DBLE(12.0))+DBLE(0.02))**(one/DBLE(12.0))) - mu ! Inflation rate ~ +4%   ---> 2% annual more
  IF (mucase.EQ."4") mu0 = (((mu**DBLE(12.0))+DBLE(0.06))**(one/DBLE(12.0))) - mu ! Inflation rate ~ +8%   ---> 6% annual more
  IF (mucase.EQ."5") mu0 = (((mu**DBLE(12.0))-DBLE(0.04))**(one/DBLE(12.0))) - mu ! Inflation rate ~ -2%   ---> 4% annual less
  IF (mucase.EQ."6") mu0 = (((mu**DBLE(12.0))-DBLE(0.01))**(one/DBLE(12.0))) - mu ! Inflation rate ~ +1%   ---> 1% annual less
  IF (mucase.EQ."7") mu0 = (               DBLE(1.046342)**(one/DBLE(12.0))) - mu ! Inflation rate ~ 4.63% ---> average 1980-2000
  IF (mucase.EQ."8") mu0 = (               DBLE(1.020054)**(one/DBLE(12.0))) - mu ! Inflation rate ~ 2.01% ---> average 2000-2020
  CALL SET_MATS( ) ! Re-fill the transition matrix to account for the chosen inflation rate
  RETURN
END SUBROUTINE SETINFLATION

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! EQUILIBRIUM
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! SUBROUTINE TO SOLVE THE STEADY STATE
SUBROUTINE SOLVESTEADY(name,IPP)
  IMPLICIT NONE
  CHARACTER(LEN=3)    , INTENT(IN) :: name
  INTEGER  , INTENT(IN) , OPTIONAL :: IPP
  INTEGER                          :: IPR
  IPR = 0 ; IF (PRESENT(IPP)) IPR = IPP
  CALL SET_MATS( )
  CALL READSTEADY(name)
  CALL COMPUTEGE(ZEROS1,PRICES1,IPR)
  CALL CALCSTATS( )
  CALL WRITESTEADY(name)
  RETURN
END SUBROUTINE SOLVESTEADY

! SUBROUTINE TO COMPUTE THE EQUAILIBRIUM
SUBROUTINE COMPUTEGE(Y1,X1,IPP)
  USE toolkit , ONLY : LMMIN
  IMPLICIT NONE
  REAL(rp) , INTENT(OUT)           :: Y1(3),X1(3)
  INTEGER  , INTENT(IN) , OPTIONAL :: IPP
  REAL(rp)                         :: X0(3),YY
  INTEGER                          :: ITER,IND,IPR
  IPR = 0
  IF (PRESENT(IPP) .AND. IPP.EQ.1) IPR = 1
  IF (PRESENT(IPP) .AND. IPP.EQ.2) IPR = 2
  IF (IPR.GT.1) WRITE(*,'(A)',ADVANCE="NO") '    Solving the steady state with initial guess... '
  Y1 = FIRSTORDERCONDITIONS(cbar,nbar,wbar)
  IF (IPR.GT.1) WRITE(*,'(A,F10.7)',ADVANCE="YES") ' Error = ', SQRT(SUM(Y1(:)*Y1(:)))
  IF (SQRT(SUM(Y1(:)*Y1(:))).GT.DBLE(0.000001)) THEN
    IF (IPR.GT.1) WRITE(*,'(A)') '    Error is too large. Computing the steady state '
    X0 = (/ cbar , nbar , wbar /)
    CALL LMMIN(EQUILIBRIUM,X1,Y1,ITER,IND,X0,ITERMAX=500,DAMP=one,SHCK=DBLE(0.05),IPRINT=IPR,USEBRO=0)
    Y1 = FIRSTORDERCONDITIONS(X1(1),X1(2),X1(3))
  END IF
  RETURN
END SUBROUTINE COMPUTEGE

! FUNCTION THAT FINDS THE STEADY STATE CONSUMPTION, WAGE AND LABOR SUPPLY
FUNCTION EQUILIBRIUM(EQPRICE) RESULT(FOCS)
  IMPLICIT NONE
  DOUBLE PRECISION              :: EQPRICE(:)
  DOUBLE PRECISION, ALLOCATABLE :: FOCS(:)
  ALLOCATE(FOCS(3))
  FOCS = FIRSTORDERCONDITIONS(EQPRICE(1),EQPRICE(2),EQPRICE(3))
  RETURN
END FUNCTION EQUILIBRIUM

! RESIDUALS FROM STEADY-STATE SYSTEM
FUNCTION FIRSTORDERCONDITIONS(cbar0,nbar0,wbar0) RESULT(FOCNDS)
  USE parameters , ONLY : nump,nums,numw,numz,Wdist,Pdist,epsilon,epsilonN,&
                          Pi,w_grid,z_grid,p_grid,s_grid,lambda
  USE solution   , ONLY : SOLVEFIRMS,SOLVEWORKERS
  IMPLICIT NONE
  DOUBLE PRECISION :: cbar0,nbar0,wbar0
  DOUBLE PRECISION :: FOCNDS(3)
  REAL(rp) :: test_c,test_n,test_w,DDelta,Kpi,Klambda,KL_pi
  INTEGER  :: is,ip,iw,iz,ips
  cbar = cbar0
  nbar = nbar0
  wbar = wbar0
  CALL SOLVEFIRMS( )
  CALL SOLVEWORKERS( )
  test_c = zero ; DDelta  = zero
  test_n = zero ; Kpi     = zero
  test_w = zero ; Klambda = zero
  ! FOC for aggregate consumption
  DO ip = 1,nump ; DO is = 1,nums
    test_c = test_c + Pdist(ip,is)*exp(p_grid(ip)*(one-epsilon))
  END DO ; END DO
  ! FOC for aggregate wgae
  DO iw = 1,numw ; DO iz = 1,numz
    test_w = test_w + Wdist(iw,iz)*exp((one-epsilonN)*(w_grid(iw)-z_grid(iz)))
  END DO ; END DO
  ! FOC for aggregate labor supply
  DO ip = 1,nump ; DO is = 1,nums
    KL_pi  = zero
    DO ips = 1,nump
      KL_pi = KL_pi + Pi(ips,is)*log(max(Pi(ips,is)*DBLE(nump),tol))
    END DO
    Kpi     = Kpi     + lambda(ip,is)*Pdist(ip,is)*KL_pi
    DDelta  = DDelta  + Pdist(ip,is)*exp(s_grid(is)-epsilon*p_grid(ip))
    Klambda = Klambda + Pdist(ip,is)*(lambda(ip,is)*log(max(tol,lambda(ip,is)/lbar)) + &
                        (one-lambda(ip,is))*log(max(tol,(one-lambda(ip,is))/(one-lbar))))
  END DO ; END DO
  test_n = DDelta*cbar + kappa_pi*Kpi + kappa_lambda*Klambda
  ! Residual in aggregate price (normalized to 1 - model in real terms)
  FOCNDS(1) = one  - test_c**(one/(one-epsilon))
  ! Residual in aggregate labor demand
  FOCNDS(2) = nbar - test_n
  ! Residual in aggregate (real) wage
  FOCNDS(3) = wbar - test_w**(one/(one-epsilonN))
  RETURN
END FUNCTION FIRSTORDERCONDITIONS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CALIBRATION
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! SUBROUTINE TO CALIBRATE THE MODEL
SUBROUTINE CALIBRATE( )
  USE toolkit    , ONLY : SIMPLEX,LMMIN,NORMALIZE
  IMPLICIT NONE
  REAL(rp) :: PARS0(7),PARS1(7),PPR(3),PPY(3)
  REAL(rp) :: TEST,TOLERANCE,SHOCK,DAMP,MOMS(144)
  INTEGER  :: ITER,ITERMAX,IPRINT,METHOD,IND

  ! Set optimization parameters
  ITERMAX   = 300                ! Maximum number of function evaluations
  TOLERANCE = DBLE(0.00000001)   ! Convergence criterium for sum of square errors.
  IPRINT    = 2                  ! 2 to print every iteration, 1 to print only important workings, 0 to not print anything
  DAMP      = 1.0000             ! Initial damping factor for Levenberg–Marquardt
  SHOCK     = DBLE(0.01)         ! Shock to compute numerical Jacobian (for Levenberg–Marquardt)

  ! Ask which method to use
  1 WRITE(*,'(A)',ADVANCE="NO") '  Calibration method: (1) Simplex, (2) LMIN : '
  READ (*,*)  METHOD
  IF ( METHOD.LT.1 .or. METHOD.GT.2) THEN
    WRITE(*,'(A)',ADVANCE="NO") '  Invalid entry...'
    GOTO 1
  END IF

  CALL READSTEADY("V10")

  ! Define moments weights in calibration
  WEIGHT = zero
  DO J=1,142
    WEIGHT(J) = one
  END DO
  WEIGHT(143) = sqrt(142.0)
  WEIGHT(144) = sqrt(142.0)

  ! Define max and min of parameters to calibrate
  PARSMAX(1) = DBLE(0.30) ; PARSMIN(1) = DBLE(0.05)    ! lbar
  PARSMAX(2) = DBLE(0.30) ; PARSMIN(2) = DBLE(0.05)    ! rhobar
  PARSMAX(3) = DBLE(0.05) ; PARSMIN(3) = DBLE(0.005)   ! kappa_w
  PARSMAX(4) = DBLE(0.97) ; PARSMIN(4) = DBLE(0.50)    ! rho_z
  PARSMAX(5) = DBLE(0.08) ; PARSMIN(5) = DBLE(0.02)    ! stdMC_z
  PARSMAX(6) = DBLE(0.97) ; PARSMIN(6) = DBLE(0.50)    ! rho_s
  PARSMAX(7) = DBLE(0.12) ; PARSMIN(7) = DBLE(0.02)    ! stdMC_s

  ! Redefine parameters to be unbounded (NORMALIZE transform bounded variables into unbounded ones, see toolkit.f90)
  CALL NORMALIZE(PARS0(1), lbar,    PARSMAX(1),PARSMIN(1),0)
  CALL NORMALIZE(PARS0(2), rhobar,  PARSMAX(2),PARSMIN(2),0)
  CALL NORMALIZE(PARS0(3), kappa_w, PARSMAX(3),PARSMIN(3),0)
  CALL NORMALIZE(PARS0(4), rho_z,   PARSMAX(4),PARSMIN(4),0)
  CALL NORMALIZE(PARS0(5), stdMC_z, PARSMAX(5),PARSMIN(5),0)
  CALL NORMALIZE(PARS0(6), rho_s,   PARSMAX(6),PARSMIN(6),0)
  CALL NORMALIZE(PARS0(7), stdMC_s, PARSMAX(7),PARSMIN(7),0)

  ! Call chosen optimization routine
  IF ( METHOD.eq.1 ) CALL SIMPLEX(SUMERRORS,PARS1,TEST,ITER,IND,PARS0,ITERMAX=1000,IPRINT=2)
  IF ( METHOD.EQ.2 ) CALL LMMIN(ERRORS,PARS1,MOMS,ITER,IND,PARS0,ITERMAX=1000,IPRINT=2)

  ! Redefine parameter values (revert normalization with NORMALIZE)
  CALL NORMALIZE(PARS0(1), lbar,    PARSMAX(1),PARSMIN(1),1)
  CALL NORMALIZE(PARS0(2), rhobar,  PARSMAX(2),PARSMIN(2),1)
  CALL NORMALIZE(PARS0(3), kappa_w, PARSMAX(3),PARSMIN(3),1)
  CALL NORMALIZE(PARS0(4), rho_z,   PARSMAX(4),PARSMIN(4),1)
  CALL NORMALIZE(PARS0(5), stdMC_z, PARSMAX(5),PARSMIN(5),1)
  CALL NORMALIZE(PARS0(6), rho_s,   PARSMAX(6),PARSMIN(6),1)
  CALL NORMALIZE(PARS0(7), stdMC_s, PARSMAX(7),PARSMIN(7),1)

  kappa_lambda = kappa_pi
  kappa_rho    = kappa_w

  CALL COMPUTEGE(PPY,PPR,2)
  CALL CALCSTATS( )
  CALL WRITESTEADY("VC0")

  RETURN
END SUBROUTINE CALIBRATE

! COMPUTE VECTOR OF CALIBRATION RESIDUALS
FUNCTION ERRORS(PARS)  RESULT(MOMS0)
  USE toolkit    , ONLY : NORMALIZE
  USE parameters , ONLY : nump,numw
  IMPLICIT NONE
  REAL(rp)               :: PARS(:),RR,PP1(3),PR1(3),vect(11)
  REAL(rp) , ALLOCATABLE :: MOMS0(:)
  INTEGER                :: j

  ! Define parameter values (revert normalization with NORMALIZE)
  CALL NORMALIZE(PARS(1), lbar,    PARSMAX(1),PARSMIN(1),1)
  CALL NORMALIZE(PARS(2), rhobar,  PARSMAX(2),PARSMIN(2),1)
  CALL NORMALIZE(PARS(3), kappa_w, PARSMAX(3),PARSMIN(3),1)
  CALL NORMALIZE(PARS(4), rho_z,   PARSMAX(4),PARSMIN(4),1)
  CALL NORMALIZE(PARS(5), stdMC_z, PARSMAX(5),PARSMIN(5),1)
  CALL NORMALIZE(PARS(6), rho_s,   PARSMAX(6),PARSMIN(6),1)
  CALL NORMALIZE(PARS(7), stdMC_s, PARSMAX(7),PARSMIN(7),1)

  kappa_lambda = kappa_pi
  kappa_rho    = kappa_w

  ! Solve for the steady-state
  CALL SET_MATS( )
  CALL COMPUTEGE(PR1,PP1,0)
  CALL CALCSTATS( )

  ! Define model residuals
  ALLOCATE(MOMS0(2*nump+2*numw))
  MOMS0(:) = DMOMS(:)

  ! Sum of squared errors
  RR = SUM(DMOMS(:)*DMOMS(:))

  ! Read last iteration's sum of squared errors
  OPEN(unit=1,file="calibparams.txt",action='read')
    DO j=1,11
      READ(1,*)  vect(j)
    END DO
  CLOSE(1)

  ! If RR is lower than the lowest found so far, write parameters
  IF (RR.lt.vect(11)) THEN
    OPEN(unit=1,file="calibparams.txt",action='write')
      WRITE(1,'(F20.15)')  lbar
      WRITE(1,'(F20.15)')  rhobar
      WRITE(1,'(F20.15)')  kappa_pi
      WRITE(1,'(F20.15)')  kappa_lambda
      WRITE(1,'(F20.15)')  kappa_w
      WRITE(1,'(F20.15)')  kappa_rho
      WRITE(1,'(F20.15)')  rho_z
      WRITE(1,'(F20.15)')  stdMC_z
      WRITE(1,'(F20.15)')  rho_s
      WRITE(1,'(F20.15)')  stdMC_s
      WRITE(1,'(F20.15)')  RR
    CLOSE (1)
    CALL CALCSTATS( )
    CALL WRITESTEADY("VC0")
  END IF

  RETURN
END FUNCTION ERRORS

! COMPUTE SUM OF SQUARED RESIDUALS FROM CALIBRATION
FUNCTION SUMERRORS(PARS)  RESULT(RR)
  IMPLICIT NONE
  REAL(rp) :: PARS(:),RR,MOMENTS(144)
  MOMENTS = ERRORS(PARS)
  RR = SUM(MOMENTS(:)*MOMENTS(:))
  RETURN
END FUNCTION SUMERRORS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! DYNAMICS
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! SUBROUTINE TO COMPUTE THE JACOBIAN OF THE DYNAMIC SUSTEM
SUBROUTINE SOLVEDYN(name)
  USE parameters , ONLY : Pdist,Wdist,V,L,nu,gamma,delta,mu,mu0,nump,nums,numz,numw
  USE dynamics   , ONLY : DYNSYS
  IMPLICIT NONE

  CHARACTER(LEN=3) , INTENT(IN) :: name
  INTEGER                       :: i,j,ip,is
  INTEGER , PARAMETER           :: NUPV = nump*nums
  INTEGER , PARAMETER           :: NUWL = numw*numz
  INTEGER , PARAMETER           :: NUWZ = 5
  INTEGER , PARAMETER           :: NUMB = 3
  INTEGER , PARAMETER           :: NUMA = 2*NUPV + 2*NUWL + NUWZ
  INTEGER , PARAMETER           :: NUMV = NUMA + NUMB
  INTEGER , PARAMETER           :: NUMR = NUMA
  REAL(rp)                      :: XVEC(NUMV*2),XVEC0(NUMV*2),PRICESD(3),ZEROSD(3),jacstep
  REAL(rp)                      :: RESIDJ(NUMR),RESID1(NUMR),RESID0(NUMR)
  REAL(rp)                      :: JACP(NUMR,NUMV),JACN(NUMR,NUMV),mbar

  jacstep = tol

  ! ----------------------------------------------------------------------------
  ! SOLVING STEADY STATE

  WRITE(*,'(A)',ADVANCE="YES") '  Solving the steady state equilibrium... '
  WRITE(*,'(A)',ADVANCE="YES") '   '

  CALL SOLVESTEADY(name,2)

  ! ----------------------------------------------------------------------------
  ! CONSTRUCT VARIABLES OF DYNAMIC SYSTEM
  ! [ Pdist', Wdist', mlag', w', Inf', V', L', C', N', ...
  !   Pdist , Wdist , mlag , w , Inf , V , L , C , N , Shocks', Shocks ]

  mbar = nu*(mu+mu0)*(cbar**gamma)/((mu+mu0)-one+delta)
  j = 0
  DO ip=1,nump ; DO is=1,nums ; j = j + 1
    XVEC0(j) = Pdist(ip,is) ; XVEC0(NUMA+j) = XVEC0(j)
  END DO ; END DO
  DO iw=1,numw ; DO iz=1,numz ; j = j + 1
    XVEC0(j) = Wdist(iw,iz) ; XVEC0(NUMA+j) = XVEC0(j)
  END DO ; END DO
  j = j + 1 ; XVEC0(j) = mbar     ; XVEC0(NUMA+j) = XVEC0(j)
  j = j + 1 ; XVEC0(j) = wbar     ; XVEC0(NUMA+j) = XVEC0(j)
  j = j + 1 ; XVEC0(j) = (mu+mu0) ; XVEC0(NUMA+j) = XVEC0(j)
  DO ip=1,nump ; DO is=1,nums ; j = j + 1
    XVEC0(j) = V(ip,is) ; XVEC0(NUMA+j) = XVEC0(j)
  END DO ; END DO
  DO iw=1,numw ; DO iz=1,numz ; j = j + 1
    XVEC0(j) = L(iw,iz) ; XVEC0(NUMA+j) = XVEC0(j)
  END DO ; END DO
  j = j + 1 ; XVEC0(j) = cbar ; XVEC0(NUMA+j) = XVEC0(j)
  j = j + 1 ; XVEC0(j) = nbar ; XVEC0(NUMA+j) = XVEC0(j)
  DO j=NUMA*2+1,NUMA*2+NUMB*2
    XVEC0(j) = zero
  END DO

  ! ----------------------------------------------------------------------------
  ! CHECKING RESIDUALS FROM STEADY STATE (SHOULD BE CLOSE TO 0)

  WRITE(*,'(A)',ADVANCE="YES") '   '
  WRITE(*,'(A)',ADVANCE="NO") '  Checking residuals in steady state... ' ; RESID0 = DYNSYS(XVEC0)
  PRINT * , MAXVAL(ABS(RESID0(:)))
  WRITE(*,'(A)',ADVANCE="YES") '   '
  WRITE(*,'(A)',ADVANCE="NO") '    Residuals in steady state: Pdist '
  PRINT * , MAXVAL(ABS(RESID0(1:NUPV)))
  WRITE(*,'(A)',ADVANCE="NO") '    Residuals in steady state: Wdist '
  PRINT * , MAXVAL(ABS(RESID0(NUPV+1:NUPV+NUWL)))
  WRITE(*,'(A)',ADVANCE="NO") '    Residuals in steady state: V     '
  PRINT * , MAXVAL(ABS(RESID0(NUPV+NUWL+2:NUPV+NUWL+NUPV+1)))
  WRITE(*,'(A)',ADVANCE="NO") '    Residuals in steady state: L     '
  PRINT * , MAXVAL(ABS(RESID0(NUPV+NUWL+NUPV+2:NUPV+NUWL+NUPV+NUWL+1)))
  WRITE(*,'(A)',ADVANCE="NO") '    Residuals in steady state: money '
  PRINT * , ABS(RESID0(NUPV+NUWL+1))
  WRITE(*,'(A)',ADVANCE="NO") '    Residuals in steady state: euler '
  PRINT * , ABS(RESID0(NUPV+NUWL+NUPV+NUWL+2))
  WRITE(*,'(A)',ADVANCE="NO") '    Residuals in steady state: labor '
  PRINT * , ABS(RESID0(NUPV+NUWL+NUPV+NUWL+3))
  WRITE(*,'(A)',ADVANCE="NO") '    Residuals in steady state: wage  '
  PRINT * , ABS(RESID0(NUPV+NUWL+NUPV+NUWL+4))
  WRITE(*,'(A)',ADVANCE="NO") '    Residuals in steady state: price '
  PRINT * , ABS(RESID0(NUPV+NUWL+NUPV+NUWL+5))
  WRITE(*,'(A)',ADVANCE="YES") '   '

  ! ----------------------------------------------------------------------------
  ! COMPUTING JACOBIAN

  WRITE(*,'(A)',ADVANCE="YES") ' Computing the jacobian... '
  WRITE(*,'(A)',ADVANCE="YES") '   '
  WRITE(*,'(A)',ADVANCE="NO" ) '    Shocking today´s variables... '
  OPEN(unit=1,file="_dyn/"//TRIM(ADJUSTL(name))//"_dyn.txt",action='write')
  DO j=1,NUMV*2
    XVEC    = XVEC0
    XVEC(j) = XVEC0(j) + jacstep*MAX(one,ABS(XVEC(j)))
    RESID1  = DYNSYS(XVEC)
    RESIDJ  = (RESID1 - RESID0)/(XVEC(j)-XVEC0(j))
    WRITE(1,*) RESIDJ(:)
    IF (j.eq.NUMV) THEN
      WRITE(*,'(A)',ADVANCE="YES") ' Done '
      WRITE(*,'(A)',ADVANCE="NO") '    Shocking tomorrow´s variables... '
    END IF
  END DO
  WRITE(*,'(A)',ADVANCE="YES") ' Done '
  CLOSE(1)

  RETURN
END SUBROUTINE SOLVEDYN

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END PROGRAM main
