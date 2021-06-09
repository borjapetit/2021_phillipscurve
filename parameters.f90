
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE parameters

IMPLICIT NONE

! Path to files  ----------------------------------------------------------------
CHARACTER(LEN=100), PARAMETER :: path = "/Users/borjapetit/Desktop/MS20191456/"
CHARACTER(LEN=13)             :: folder

! Program execution parameters -------------------------------------------------
INTEGER  , PARAMETER :: rp     = kind(1.0d00)   ! Precision for real numbers.
REAL(rp) , PARAMETER :: zero   = DBLE(0.0000)   ! 0 in double precision.
REAL(rp) , PARAMETER :: one    = DBLE(1.0000)   ! 1 in double precision.
REAL(rp) , PARAMETER :: diez   = DBLE(10.000)   ! 10 in double precision.
REAL(rp) , PARAMETER :: cien   = DBLE(100.00)   ! 100 in double precision.
REAL(rp) , PARAMETER :: mil    = DBLE(1000.0)   ! 1000 in double precision.
INTEGER              :: numthreads              ! Dummy for the number of threads when using OpenMP

! Tolerance level --------------------------------------------------------------
REAL(rp) , PARAMETER :: tol    = 0.014901161193847656/(mil*mil) ! Tol value in Matlab

! Aggregate shocks and policy parameters ---------------------------------------
REAL(rp) , PARAMETER :: mu    = (one+4.2604e-004)**DBLE(4.0)  ! Monthly steady state rate of money growth (roughly 2%/yr)
REAL(rp) , PARAMETER :: rho_r = 0.8000                        ! Monthly persistence of money shock
REAL(rp)             :: mu0                                   ! Monthly steady state rate of money growth (extra - for experiments)

! Labor Productivity process ---------------------------------------------------
REAL(rp) , PARAMETER :: SSProd    =  0.3         ! Average workers' productivity for Tauchen's discretization
REAL(rp) , PARAMETER :: InitProd  = -0.3         ! Initial workers' productivity
REAL(rp) , PARAMETER :: deathprob = 0.025/12.0   ! Monthly workers' death probability (lifetime ~40yr)

! Discount and preference parameters -------------------------------------------
REAL(rp) , PARAMETER :: delta     = one-DBLE(1.04)**(-one/DBLE(12.0))  ! MONTHLY time discount factor
REAL(rp) , PARAMETER :: gamma     = 2.00000                            ! CRRA exponent
REAL(rp) , PARAMETER :: chi       = 6.00000                            ! Labor supply coefficient
REAL(rp) , PARAMETER :: labelast  = 0.50000                            ! Labor supply exponent: labelast=0 implies LINEAR labor disutility
REAL(rp) , PARAMETER :: nu        = 1.00000                            ! Money demand parameter
REAL(rp) , PARAMETER :: epsilon   = 7.00000                            ! Elasticity of substitution among product varieties
REAL(rp) , PARAMETER :: epsilonN  = 7.00000                            ! Elasticity of substitution among workers

! State Space dimension --------------------------------------------------------
INTEGER , PARAMETER :: nump = 36   ! Number of points in the grid for prices
INTEGER , PARAMETER :: nums = 25   ! Number of points in the grid for firms' productivity
INTEGER , PARAMETER :: numw = 71   ! Number of points in the grid for wages
INTEGER , PARAMETER :: numz = 51   ! Number of points in the grid for workers' productivity

! State variables --------------------------------------------------------------
REAL(rp) :: p_grid(nump)    ! Grid for (log) prices
REAL(rp) :: w_grid(numw)    ! Grid for (log) wages
REAL(rp) :: s_grid(nums)    ! Grid for (log) firms' productivity
REAL(rp) :: z_grid(numz)    ! Grid for (log) workers' productivity

! Transition matrices ----------------------------------------------------------
REAL(rp) :: T(nump,nump)      ! Transition matrix for real prices
REAL(rp) :: Tw(numw,numw)     ! Transition matrix for real wages
REAL(rp) :: S_s(nums,nums)    ! Transition matrix for firms' productivity
REAL(rp) :: S_z(numz,numz)    ! Transition matrix for workers' productivity

! Firms problem ----------------------------------------------------------------
REAL(rp) :: V(nump,nums)          ! Value function
REAL(rp) :: Pdist(nump,nums)      ! Distribution
REAL(rp) :: Pi(nump,nums)         ! Next period's price distribution
REAL(rp) :: lambda(nump,nums)     ! Probability of adjutment
REAL(rp) :: LMU(nump,nums)        ! TIme devoted to timing choice
REAL(rp) :: LABUS(nump,nums)      ! Total decision time (timing + size)

! Workers problem --------------------------------------------------------------
REAL(rp) :: L(numw,numz)          ! Value function
REAL(rp) :: Wdist(numw,numz)      ! Distribution
REAL(rp) :: Pi_w(numw,numz,numw)  ! Next period's wage distribution
REAL(rp) :: rho(numw,numz)        ! Probability of adjustment

REAL(rp) :: BirthDist(numw,numz)
REAL(rp) :: Xp(numw,numz)
REAL(rp) :: H(numw,numz)
REAL(rp) :: Ld(numw,numz)
REAL(rp) :: Lrho(numw,numz)
REAL(rp) :: AUXMAT(numw,numz)

! Equilibrium "prices" ---------------------------------------------------------
REAL(rp) :: wbar      ! Steady state wage
REAL(rp) :: cbar      ! Steady state consumption
REAL(rp) :: nbar      ! Steady state labor supply

! Cost parameters --------------------------------------------------------------
REAL(rp) :: lbar
REAL(rp) :: rhobar
REAL(rp) :: kappa_pi
REAL(rp) :: kappa_lambda
REAL(rp) :: kappa_w
REAL(rp) :: kappa_rho

! Labor productivity parameters (to be estimated)  -----------------------------
REAL(rp) :: stdMC_z              ! Standard deviation of workers' productivity process
REAL(rp) :: rho_z                ! MONTHLY persistence of workers' productivity shock process
REAL(rp) :: stdMC_s              ! Standard deviation of firms' productivity process
REAL(rp) :: rho_s                ! MONTHLY persistence of firms' productivity shock process

! Empirical moments (Data & Model)  --------------------------------------------
REAL(rp) :: DMOMS(2*nump+2*numw)
REAL(rp) :: MOMMODEL(2*nump+2*numw)
REAL(rp) :: MOMDATA(2*nump+2*numw)
REAL(rp) :: WEIGHT(2*nump+2*numw)
REAL(rp) :: MOMPRINTMODEL(38)
REAL(rp) :: MOMPRINTDATA(38)
CHARACTER(LEN=15) :: MOMPRINTNAME(19)

CONTAINS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE SET_MATS( )
    USE toolkit  , ONLY : GRID,TAUCHEN,TRANSMAT
    IMPLICIT NONE
    INTEGER  :: j
    REAL(rp) :: Smat(nums,nums),Zmat(numz,numz)

    ! Define grids for state variables.
    p_grid = GRID(DBLE(0.35),-DBLE(0.35),nump,one)
    w_grid = GRID(DBLE(0.35),-DBLE(0.35),numw,one)
    s_grid = GRID(DBLE(4.0)*stdMC_s,-DBLE(4.0)*stdMC_s,nums,one)
    z_grid = GRID(DBLE(4.0)*stdMC_z+SSProd,-DBLE(2.0)*stdMC_z+InitProd,numz,one)

    ! Define matrix T in steady-state
    T  = RMATRIX(mu+mu0)
    Tw = RMATRIXw(mu+mu0)

    ! Define transition matrices for shocks s and z
    CALL TAUCHEN(s_grid,rho_s,zero                  ,stdMC_s*SQRT(one-rho_s*rho_s),nums,Smat) ; S_s = TRANSMAT(Smat)
    CALL TAUCHEN(z_grid,rho_z,SSProd*DBLE(one-rho_z),stdMC_z*SQRT(one-rho_z*rho_z),numz,Zmat) ; S_z = TRANSMAT(Zmat)

    ! Define moments weights in calibration
    WEIGHT = zero
    DO J=1,142
      WEIGHT(J) = one
    END DO
    WEIGHT(143) = sqrt(142.0)
    WEIGHT(144) = sqrt(142.0)

    RETURN
  END SUBROUTINE SET_MATS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  FUNCTION RMATRIX(INF) RESULT(TMAT)
    IMPLICIT NONE
    REAL(rp) :: INF,TMAT(nump,nump),aux,aux1
    INTEGER  :: i
    TMAT = zero
    aux  = LOG(INF)/(DBLE(0.35)*DBLE(2.0)/(nump-one))
    aux1 = aux - DBLE(FLOOR(aux))
    IF (aux.GE.zero) THEN
      DO i = 1,nump-1
        TMAT(i,i)   = one - aux1
        TMAT(i,i+1) = aux1
      END DO
      TMAT(nump,nump) = one-aux
      TMAT(1,1)       = one
    ELSE
      DO i = 1,nump-1
        TMAT(i,i)   = aux1
        TMAT(i+1,i) = one-aux1
      END DO
      TMAT(nump,nump) = one
    END IF
    RETURN
  END FUNCTION RMATRIX
  FUNCTION RMATRIXw(INF) RESULT(TMAT)
    IMPLICIT NONE
    REAL(rp) :: INF,TMAT(numw,numw),aux,aux1
    INTEGER  :: i
    TMAT = zero
    aux  = LOG(INF)/(DBLE(0.35)*DBLE(2.0)/(numw-one))
    aux1 = aux - DBLE(FLOOR(aux))
    IF (aux.GE.zero) THEN
      DO i = 1,numw-1
        TMAT(i,i)   = one - aux1
        TMAT(i,i+1) = aux1
      END DO
      TMAT(numw,numw) = one-aux
      TMAT(1,1)       = one
    ELSE
      DO i = 1,numw-1
        TMAT(i,i)   = aux1
        TMAT(i+1,i) = one-aux1
      END DO
      TMAT(numw,numw) = one
    END IF
    RETURN
  END FUNCTION RMATRIXw

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  FUNCTION BIRTHDISTRIBUTION(consumption,wagerate,agglabor) RESULT(BDIST)
    USE toolkit , ONLY : INTERPOLATION
    IMPLICIT NONE
    REAL(rp) :: BDIST(numw,numz),consumption,wagerate,agglabor,D0(numw),relt_w,flexw_iz,Dz(numz)
    INTEGER  :: iw,iz
    CALL INTERPOLATION(iw,relt_w,InitProd,z_grid)
    D0    = zero ; D0(iw) = relt_w ; D0(iw-1) = one-relt_w
    BDIST = zero ; Dz     = zero
    DO iz=1,numz
      flexw_iz = ( log(epsilonN*chi/(epsilonN-one)) + gamma*log(consumption) + &
                 labelast*log(agglabor) + (labelast*epsilonN)*log(wagerate) )/(one+labelast*epsilonN) + &
                 (labelast*(epsilonN-one)/(one+labelast*epsilonN))*z_grid(iz)
      Dz(iz) = flexw_iz
      CALL INTERPOLATION(iw,relt_w,flexw_iz,w_grid)
      BDIST(iw,iz)   = relt_w*D0(iz)
      BDIST(iw-1,iz) = (one-relt_w)*D0(iz)
    END DO
    RETURN
  END FUNCTION BIRTHDISTRIBUTION

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE SET_MOMENTS
    IMPLICIT NONE
    INTEGER  :: i,j

    OPEN(unit=1,file="data_pdfprices.txt",action='read')
      DO i = 1, (2*nump-1) ; READ(1,*)  MOMDATA(i) ; END DO
    CLOSE (1)
    OPEN(unit=1,file="data_pdfwages.txt",action='read') ; j = 0
      DO i = (2*nump),(2*nump+2*numw-2) ; j = j + 1 ; READ(1,*)  MOMDATA(i) ; END DO
    CLOSE (1)
    MOMDATA(2*nump+2*numw-1) = 0.102
    MOMDATA(2*nump+2*numw)   = 1.0/12.0

    MOMPRINTDATA(1)  = 10.2  ; MOMPRINTDATA(20) = 8.33
    MOMPRINTDATA(2)  = 1.60  ; MOMPRINTDATA(21) = 5.10
    MOMPRINTDATA(3)  = 9.90  ; MOMPRINTDATA(22) = 6.47
    MOMPRINTDATA(4)  = 13.2  ; MOMPRINTDATA(23) = 6.52
    MOMPRINTDATA(5)  = -0.42 ; MOMPRINTDATA(24) = 0.35
    MOMPRINTDATA(6)  = 4.81  ; MOMPRINTDATA(25) = 4.39
    MOMPRINTDATA(7)  = 65.1  ; MOMPRINTDATA(26) = 86.5
    MOMPRINTDATA(8)  = 35.5  ; MOMPRINTDATA(27) = 43.0
    MOMPRINTDATA(9)  = 12.0  ; MOMPRINTDATA(28) = 11.8
    MOMPRINTDATA(10) = zero  ; MOMPRINTDATA(29) = zero
    MOMPRINTDATA(11) = zero  ; MOMPRINTDATA(30) = zero
    MOMPRINTDATA(12) = zero  ; MOMPRINTDATA(31) = zero
    MOMPRINTDATA(13) = zero  ; MOMPRINTDATA(32) = zero
    MOMPRINTDATA(14) = zero  ; MOMPRINTDATA(33) = zero
    MOMPRINTDATA(15) = zero  ; MOMPRINTDATA(34) = zero

    RETURN
  END SUBROUTINE SET_MOMENTS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  FUNCTION JIMDIAG(A,J) RESULT(diagvec)
    INTEGER :: L,J,K,startrow,startcol,diaglen
    REAL(rp) :: A(:,:)
    REAL(rp), ALLOCATABLE :: diagvec(:)
    K       = size(A,1)
    diaglen = K-abs(J)
    ALLOCATE(diagvec(diaglen))
    IF (J.eq.0) THEN
     startrow = 1
     startcol = 1
    ELSE IF (J.gt.0) THEN
      startrow = 1
      startcol = J
    ELSE IF (J.lt.0) THEN
      startrow = abs(J)
      startcol = 1
    END IF
    DO L = 1,diaglen
      diagvec(L) = A(startrow+L-1,startcol+L-1)
    END DO
    RETURN
  END FUNCTION JIMDIAG

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE PRINTVERSION(versionname)
    IMPLICIT NONE
    CHARACTER(LEN=3) , INTENT(IN)  :: versionname
    CHARACTER(LEN=22)              :: verskappa
    IF (versionname(2:2).EQ."1") THEN ; verskappa = " | ST p | ST w | pi = " ; END IF
    IF (versionname(2:2).EQ."2") THEN ; verskappa = " | SF p | ST w | pi = " ; END IF
    IF (versionname(2:2).EQ."3") THEN ; verskappa = " | FL p | ST w | pi = " ; END IF
    IF (versionname(2:2).EQ."4") THEN ; verskappa = " | ST p | SF w | pi = " ; END IF
    IF (versionname(2:2).EQ."5") THEN ; verskappa = " | ST p | FL w | pi = " ; END IF
    IF (versionname(2:2).EQ."6") THEN ; verskappa = " | FL p | FL w | pi = " ; END IF
    IF (versionname(2:2).NE."C") THEN
      WRITE(*,'(A)',ADVANCE="NO") '  Loading initial guess... '
      WRITE(*,'(A,I2)') ' Version '//versionname//verskappa , NINT(cien*((mu+mu0)**(12.00)) - cien)
      WRITE(*,'(A)') '  '
    END IF
    RETURN
  END SUBROUTINE PRINTVERSION

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE READSTEADY(name)
    IMPLICIT NONE
    LOGICAL :: OK,OK0
    CHARACTER(LEN=3) , INTENT(IN), OPTIONAL :: name

    ! If initial guess available, load it
    IF (PRESENT(name)) THEN
      INQUIRE(file="_ss/"//TRIM(ADJUSTL(name))//"_ss.txt",EXIST=OK)
      IF (OK) THEN
        CALL READFILE(name)
        RETURN
      END IF
    END IF

    ! If no initial guess available, load the one corresponding with inflation 2% keeping adjustment costs
    INQUIRE(file="_ss/V"//name(2:2)//"0_ss.txt",EXIST=OK0)
    IF (OK0) THEN
      CALL READFILE("V"//name(2:2)//"0")
      WRITE(*,'(A,A,A)') '    No initial guess available. Setting V',name(2:2),'0 as initial guess'
      RETURN
    END IF

    ! If no initial guess available, load the one corresponding with baseline scenario
    INQUIRE(file="_ss/V10_ss.txt",EXIST=OK0)
    IF (OK0) THEN
      CALL READFILE("V10")
      WRITE(*,'(A,A,A)') '    No initial guess available. Setting V1 as initial guess'
      RETURN
    END IF

    ! If no initial guess, set a default one
    CALL READDEFAULT( )
    WRITE(*,'(A)') '    No initial guess available. Setting default initial guess'

    RETURN
    CONTAINS
      SUBROUTINE READFILE(name2)
        IMPLICIT NONE
        INTEGER :: iw,ip,iz,is
        CHARACTER(LEN=3) , INTENT(IN) :: name2
        OPEN(unit=4,file="_ss/"//TRIM(ADJUSTL(name2))//"_ss.txt",action='READ')
          READ(4,*) cbar
          READ(4,*) nbar
          READ(4,*) wbar
          DO ip=1,nump ; DO is=1,nums ; READ(4,*) Pdist(ip,is) ; END DO ; END DO
          DO ip=1,nump ; DO is=1,nums ; READ(4,*) V(ip,is)     ; END DO ; END DO
          DO iw=1,numw ; DO iz=1,numz ; READ(4,*) Wdist(iw,iz) ; END DO ; END DO
          DO iw=1,numw ; DO iz=1,numz ; READ(4,*) L(iw,iz)     ; END DO ; END DO
        CLOSE(4)
        RETURN
      END SUBROUTINE READFILE
      SUBROUTINE READDEFAULT( )
        IMPLICIT NONE
        INTEGER :: iw,ip,iz,is
        cbar  = 0.50
        nbar  = 0.50
        wbar  = 0.85
        Pdist = one/DBLE(nump*nums)
        Wdist = one/DBLE(numw*numz)
        DO ip = 1,nump ; DO is = 1,nums
          V(ip,is) = (one/delta)*( cbar*exp(p_grid(ip)*(one-epsilon)) - wbar*cbar*exp(-p_grid(ip)*epsilon)*exp(s_grid(is)) )
        END DO ; END DO
        DO iw = 1,numw ; DO iz = 1,numz
          L(iw,iz) = (one/delta)*( exp(w_grid(iw))*(wbar**epsilonN)*nbar*exp( z_grid(iz)*(epsilonN-one) - w_grid(iw)*epsilonN ) &
                      - chi*(Ld(iw,iz)**(one+labelast))/((cbar**(-gamma))*(one+labelast)) )
        END DO ; END DO
        RETURN
      END SUBROUTINE READDEFAULT
  END SUBROUTINE READSTEADY

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE WRITESTEADY(name0)
    IMPLICIT NONE
    INTEGER :: iw,ip,iz,is,iwe
    CHARACTER(LEN=3) , INTENT(IN), OPTIONAL :: name0
    CHARACTER(LEN=3) :: name
    name = "V00" ; IF (PRESENT(name0)) name = name0(:)
    OPEN(unit=4,file="_ss/"//TRIM(ADJUSTL(name))//"_ss.txt",action='write')
      WRITE(4,*) cbar
      WRITE(4,*) nbar
      WRITE(4,*) wbar
      DO ip=1,nump ; DO is=1,nums ; WRITE(4,*) Pdist(ip,is)      ; END DO ; END DO
      DO ip=1,nump ; DO is=1,nums ; WRITE(4,*) V(ip,is)          ; END DO ; END DO
      DO iw=1,numw ; DO iz=1,numz ; WRITE(4,*) Wdist(iw,iz)      ; END DO ; END DO
      DO iw=1,numw ; DO iz=1,numz ; WRITE(4,*) L(iw,iz)          ; END DO ; END DO
      DO iw=1,SIZE(MOMMODEL)      ; WRITE(4,*) MOMMODEL(iw)      ; END DO
      DO iw=1,SIZE(MOMPRINTMODEL) ; WRITE(4,*) MOMPRINTMODEL(iw) ; END DO
    CLOSE(4)
    IF (name.EQ."V10") THEN
      OPEN(unit=4,file="_ss/_V10_lambda_ss.txt",action='write')
        DO ip=1,nump ; DO is=1,nums ; WRITE(4,*) lambda(ip,is) ; END DO ; END DO
      CLOSE(4)
      OPEN(unit=4,file="_ss/_V10_pi_ss.txt",action='write')
        DO ip=1,nump ; DO is=1,nums ; WRITE(4,*) Pi(ip,is) ; END DO ; END DO
      CLOSE(4)
      OPEN(unit=4,file="_ss/_V10_rho_ss.txt",action='write')
        DO iw=1,numw ; DO iz=1,numz ; WRITE(4,*) rho(iw,iz) ; END DO ; END DO
      CLOSE(4)
      OPEN(unit=4,file="_ss/_V10_piw_ss.txt",action='write')
        DO iw=1,numw ; DO iz=1,numz ; DO iwe=1,numw ; WRITE(4,*) Pi_w(iw,iz,iwe) ; END DO ; END DO ; END DO
      CLOSE(4)
    END IF
    RETURN
  END SUBROUTINE WRITESTEADY

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This subroutine computes a set of statistics from the model, including the
! model-generated moments used in calibration

  SUBROUTINE CALCSTATS( )
    USE toolkit , ONLY : GRID,STATS,VECTORIZE
    IMPLICIT NONE
    INTEGER  :: ip,is,iw,iz,ips,iws,j,j1,j2
    REAL(rp) :: AveP,AveW,StdP,StdW,Avelambda,Averho,Avelambdapos,Averhopos
    REAL(rp) :: AvePCh,AveWCh,StdPCh,StdWCh,SkePCh,SkeWCh
    REAL(rp) :: KurPCh,KurWCh,AvePChAbs,AveWChAbs
    REAL(rp) :: AvePstar,AveWstar,StdPstar,StdWstar
    REAL(rp) :: masspchanges(2*nump-1),masswchanges(2*numw-1),vecpchanges(2*nump-1)
    REAL(rp) :: denspchanges(2*nump-1),denswchanges(2*numw-1),vecwchanges(2*numw-1)
    REAL(rp) :: pStar(nums),wStar(numz),Pdergo(nump,nums),Wdergo(numw,numz)
    REAL(rp) :: LossinOptProf_p,LossinOptRev_p,LossinRev_p
    REAL(rp) :: LossinOptProf_w,LossinOptRev_w,LossinRev_w
    REAL(rp) :: aux1,aux2,aux3,aux4,aux5
    REAL(rp) :: atilde,atildeF,ztilde,ztildeF,zhat,relwF(numz)
    REAL(rp) :: Htot,Ntot,mutot,tautot,muwtot,tauwtot
    REAL(rp) :: totloss ,totloss_1 ,totloss_2,totloss_3,totloss_4
    REAL(rp) :: totlossp,totlossp_1,totlossp_2
    REAL(rp) :: totlossw,totlossw_1,totlossw_2
    REAL(rp) :: aggNF,aggHF

    ! Mean and STD of prices and wages (in logs)
    AveP = zero ; AveW = zero
    StdP = zero ; StdW = zero
    DO ip=1,nump
      AveP = AveP + p_grid(ip)*sum(Pdist(ip,:))     ! Average price
      AveW = AveW + w_grid(ip)*sum(Wdist(ip,:))     ! Average wage
    END DO
    DO ip=1,nump
      StdP = StdP + ((p_grid(ip) - AveP)**DBLE(2.00))*sum(Pdist(ip,:))
      StdW = StdW + ((w_grid(ip) - AveW)**DBLE(2.00))*sum(Wdist(ip,:))
    END DO
    StdP = StdP**DBLE(0.5)    ! Standard deviaition of prices
    StdW = StdW**DBLE(0.5)    ! Standard deviaition of wages

    ! Desity of price changes and frequency of adjustment, model
    vecpchanges  = GRID(p_grid(nump)-p_grid(1),p_grid(1)-p_grid(nump),2*nump-1,one)
    masspchanges = zero
    DO is=1,nums ; DO ip = 1,nump
      DO ips = 1,ip-1
      masspchanges(nump-ip+ips) = masspchanges(nump-ip+ips) + lambda(ip,is)*Pdist(ip,is)*Pi(ips,is)
      END DO
      masspchanges(nump)        = masspchanges(nump)        + lambda(ip,is)*Pdist(ip,is)*Pi(ips,is)
      DO ips = ip+1,nump
      masspchanges(nump+ips-ip) = masspchanges(nump+ips-ip) + lambda(ip,is)*Pdist(ip,is)*Pi(ips,is)
      END DO
    END DO ; END DO
    denspchanges = masspchanges/SUM(masspchanges)                         ! Histogram of price changes
    Avelambda    = SUM(masspchanges)                                      ! Average price adjustment probability
    Avelambdapos = SUM(masspchanges(nump+1:2*nump-1))/SUM(masspchanges)   ! Share of price increases

    ! Desity of wage changes and frequency of adjustment, model
    vecwchanges  = GRID(w_grid(numw)-w_grid(1),w_grid(1)-w_grid(numw),2*numw-1,one)
    masswchanges = zero
    DO iz=1,numz ; DO iw = 1,numw
      DO iws = 1,iw-1
      masswchanges(numw-iw+iws) = masswchanges(numw-iw+iws) + rho(iw,iz)*Wdist(iw,iz)*Pi_w(iw,iz,iws)
      END DO
      masswchanges(numw)        = masswchanges(numw)        + rho(iw,iz)*Wdist(iw,iz)*Pi_w(iw,iz,iws)
      DO iws = iw+1,numw
      masswchanges(numw+iws-iw) = masswchanges(numw+iws-iw) + rho(iw,iz)*Wdist(iw,iz)*Pi_w(iw,iz,iws)
      END DO
    END DO ; END DO
    denswchanges = masswchanges/SUM(masswchanges)                         ! Histogram of wage changes
    Averho       = SUM(masswchanges)                                      ! Average wage adjustment probability
    Averhopos    = SUM(masswchanges(numw+1:2*numw-1))/SUM(masswchanges)   ! Share of wage increases

    ! Model moments for calibration
    MOMMODEL(1:2*nump-1)             = denspchanges   ! Histogram of price changes
    MOMMODEL(2*nump:2*nump+2*numw-2) = denswchanges   ! Histogram of wage changes
    MOMMODEL(2*nump+2*numw-1)        = Avelambda      ! Average price adjustment probaility
    MOMMODEL(2*nump+2*numw)          = AveRho         ! Average wage adjustment probaility

    ! Moments deviation
    DMOMS = zero
    DO j=1,(2*nump+2*numw)
      DMOMS(j) = WEIGHT(j)*( MOMMODEL(j) - MOMDATA(j) )*DBLE(10.00)
    END DO
    DMOMS(2*nump+2*numw-1) = WEIGHT(2*nump+2*numw-1)*( MOMMODEL(2*nump+2*numw-1)/MOMDATA(2*nump+2*numw-1) - one )*DBLE(10.00)
    DMOMS(2*nump+2*numw)   = WEIGHT(2*nump+2*numw  )*( MOMMODEL(2*nump+2*numw  )/MOMDATA(2*nump+2*numw  ) - one )*DBLE(10.00)

    ! **************************************************************************
    ! Compute toher statiscts of price/wage changes

    ! Compute statistics of price changes
    CALL STATISTICS(AvePCh,AvePChAbs,StdPCh,SkePCh,KurPCh,vecpchanges,denspchanges)

    ! Compute statistics of wage changes
    CALL STATISTICS(AveWCh,AveWChAbs,StdWCh,SkeWCh,KurWCh,vecwchanges,denswchanges)

    ! **************************************************************************
    ! Losses relative to flexible price setting

    ! Optimal price level
    DO is=1,nums
      pStar(is) = LOG(EXP(s_grid(is))*epsilon/(epsilon-one)*wbar)
    END DO

    ! Average and STD of price deviations from optimal
    AvePstar = zero ; StdPstar = zero
    DO is=1,nums ; DO ip = 1,nump
      AvePstar = AvePstar + abs(p_grid(ip)-pStar(is))*Pdist(ip,is)
    END DO ; END DO
    DO is=1,nums ; DO ip = 1,nump
      StdPstar = StdPstar + ((abs(p_grid(ip)-pStar(is)) - AvePstar)**DBLE(2.00))*Pdist(ip,is)
    END DO ; END DO
    StdPstar = sqrt(StdPstar)

    ! Ergodic distribution of prices
    DO is=1,nums ; Pdergo(:,is) = zero
      j1 = 0 ; j2 = 0
      DO WHILE (j2.EQ.0) ; j1 = j1 + 1
        IF (p_grid(j1)>pStar(is)) j2 = j1
        IF (j1.eq.nump          ) j2 = nump
      END DO
      j2 = max(1,j1-1)
      Pdergo(j1,is) = SUM(Pdist(:,is))*(pStar(is)-p_grid(j2))/(p_grid(2)-p_grid(1))
      Pdergo(j2,is) = SUM(Pdist(:,is))*(p_grid(j1)-pStar(is))/(p_grid(2)-p_grid(1))
    END DO

    aux1 = zero ; aux2 = zero ; aux3 = zero ; aux4 = zero ; aux5 = zero
    DO ip = 1,nump ; DO is=1,nums
      aux1 = aux1 + Pdist(ip,is)*cbar*exp(p_grid(ip)*(one-epsilon))
      aux2 = aux2 + Pdist(ip,is)*wbar*cbar*exp(-p_grid(ip)*epsilon)*exp(s_grid(is))
      aux3 = aux3 + Pdist(ip,is)*wbar*LABUS(ip,is)
      aux4 = aux4 + Pdergo(ip,is)*cbar*exp(p_grid(ip)*(one-epsilon))
      aux5 = aux5 + Pdergo(ip,is)*wbar*cbar*exp(-p_grid(ip)*epsilon)*exp(s_grid(is))
    END DO ; END DO
    LossinOptProf_p = ( (aux4-aux5) - (aux1-aux2-aux3) ) / (aux4-aux5)    ! Losses as % optimal profits
    LossinOptRev_p  = ( (aux4-aux5) - (aux1-aux2-aux3) ) / aux4           ! Losses as % optimal revenues
    LossinRev_p     = ( (aux4-aux5) - (aux1-aux2-aux3) ) / aux1           ! Losses as % revenues

    ! **************************************************************************
    ! Losses relative to flexible wage setting

    ! Optimal wage (conditional of productivity)
    DO iz=1,numz
      wStar(iz) = (one/(one+epsilonN*labelast))*(LOG(((epsilonN*chi)/(epsilonN-one))*&
      (cbar**gamma)*(nbar**labelast)*(wbar**(epsilonN*labelast))) &
      + labelast*(epsilonN-one)*z_grid(iz))
    END DO

    ! Average and STD of wage deviations from optimal
    AveWstar = zero ; StdWstar = zero
    DO iz=1,numz ; DO iw = 1,numw
      AveWstar = AveWstar + abs(w_grid(iw)-wStar(iz))*Wdist(ip,iz)
    END DO ; END DO
    DO iz=1,numz ; DO iw = 1,numw
      StdWstar = StdWstar + ((abs(w_grid(iw)-wStar(iz)) - AveWstar)**DBLE(2.00))*Wdist(iw,iz)
    END DO ; END DO
    StdWstar = sqrt(StdWstar)

    ! Ergodic distribution of wages
    DO iz=1,numz ; Wdergo(:,iz) = zero
      j1 = 0 ; j2 = 0
      DO WHILE (j2.EQ.0) ; j1 = j1 + 1
        IF (w_grid(j1)>wStar(iz)) j2 = j1
        IF (j1.eq.numw          ) j2 = numw
      END DO
      j2 = max(1,j1-1)
      Wdergo(j1,iz) = SUM(Wdist(:,iz))*(wStar(iz)-w_grid(j2))/(w_grid(2)-w_grid(1))
      Wdergo(j2,iz) = SUM(Wdist(:,iz))*(w_grid(j1)-wStar(iz))/(w_grid(2)-w_grid(1))
    END DO

    aux1 = zero ; aux2 = zero ; aux3 = zero ; aux4 = zero ; aux5 = zero
    DO iw = 1,numw ; DO iz=1,numz
      aux1 = aux1 + Wdist(iw,iz)*exp(w_grid(iw))*Ld(iw,iz)
      aux2 = aux2 + Wdist(iw,iz)*chi*((H(iw,iz)**(one+labelast))/(one+labelast))*(cbar**gamma)
      aux3 = aux3 + Wdergo(iw,iz)*exp(w_grid(iw))*Ld(iw,iz)
      aux4 = aux4 + Wdergo(iw,iz)*chi*((Ld(iw,iz)**(one+labelast))/(one+labelast))*(cbar**gamma)
      aux5 = aux5 + Wdergo(iw,iz)*exp(w_grid(iw))*Ld(iw,iz)
    END DO ; END DO
    LossinOptProf_w = ( (aux3-aux4) - (aux1-aux2) ) / (aux3-aux4)    ! Losses as % optimal labor income net utility cost
    LossinOptRev_w  = ( (aux3-aux4) - (aux1-aux2) ) / aux5           ! Losses as % optimal labor income
    LossinRev_w     = ( (aux3-aux4) - (aux1-aux2) ) / aux1           ! Losses as % labor income

    ! **************************************************************************
    ! Cost of frictions

    Ntot    = Nbar                                                        ! Total labor input
    mutot   = SUM(Pdist(:,:)*LMU(:,:))                                    ! Labor input used for timing choice
    tautot  = SUM(Pdist(:,:)*LABUS(:,:)) - mutot                          ! Labor input used for setting choice

    Htot    = SUM(Wdist(:,:)*H(:,:))                                      ! Total labor effort
    muwtot  = SUM(Wdist(:,:)*LRHO(:,:))                                   ! Total labor effort used for timing choice
    tauwtot = SUM(Wdist(:,:)*H(:,:)) - SUM(Wdist(:,:)*Ld(:,:)) - muwtot   ! Total labor effort used for setting choice

    ! Equation 37
    atilde = zero
    DO is=1,nums ; DO ip = 1,nump
      atilde = atilde + Pdist(ip,is)*EXP(-epsilon*p_grid(ip)+s_grid(is))
    END DO ; END DO
    atilde  = one/atilde

    ! Euqation 47
    atildeF = zero
    DO is=1,nums ; DO ip = 1,nump
      atildeF = atildeF + Pdist(ip,is)*EXP(-s_grid(is)*(epsilon-one))
    END DO ; END DO
    atildeF  = atildeF**(one/(epsilon-one))

    ! equation 39
    ztilde = zero
    DO iz = 1,numz ; DO iw = 1,numw
      ztilde = ztilde + Wdist(iw,iz)*EXP(z_grid(iz)*(epsilonN-one))*((EXP(w_grid(iw))/wbar)**(-epsilonN))
    END DO ; END DO
    ztilde  = one/ztilde

    ! Equations 52 and 53
    zhat = zero
    DO iz = 1,numz ; DO iw = 1,numw
      zhat = zhat + Wdist(iw,iz)*EXP(z_grid(iz)*(one+labelast)*(epsilonN-one)/(one+labelast*epsilonN))
    END DO ; END DO
    DO iz = 1,numz
      relwF(iz) = EXP(z_grid(iz)*labelast*(epsilonN-one)/(one+labelast*epsilonN))*(zhat**(one/(epsilonN-one)))
    END DO
    ztildeF = zero
    DO iz = 1,numz ; DO iw = 1,numw
      ztildeF = ztildeF + Wdist(iw,iz)*EXP(z_grid(iz)*(epsilonN-one))*(relwF(iz)**(-epsilonN))
    END DO ; END DO
    ztildeF  = one/ztildeF

    ! Equation 79
    aggNF = ((atildeF**(one-gamma))*((epsilon-one)/epsilon)*((epsilonN-one)/epsilonN)*(one/chi)*&
            (zhat**((one+labelast*epsilonN)/(epsilonN-one))))**(one/(gamma+labelast))
    aggHF = aggNF/ztildeF

    ! **************************************************************************
    ! TOTAL LOSSES (Equation 44)

    ! Total losses (equation 44)
    totloss_1 = (atildeF*ztildeF-atilde*ztilde)/(atildeF*ztildeF)                           ! Total losses - Misallocation
    totloss_2 = (atilde*ztilde                )/(atildeF*ztildeF)*(muwtot+tauwtot)/aggHF    ! Total losses - Wage setting
    totloss_3 = (atilde                       )/(atildeF*ztildeF)*(mutot+tautot)/aggHF      ! Total losses - Price setting
    totloss_4 = (atilde*ztilde                )/(atildeF*ztildeF)*(aggHF-Htot)/aggHF        ! Total losses - Labor supply
    totloss   = totloss_1 + totloss_2 + totloss_3 + totloss_4                               ! Total losses

    ! Output losses due to price stickiness (equation 42)
    totlossp_1 = (atildeF-atilde)/(atildeF)
    totlossp_2 = (atilde)/(atildeF)*(mutot+tautot)/Ntot
    totlossp   = totlossp_1 + totlossp_2

    ! Output losses due to wage stickiness (equation 43)
    totlossw_1 = (ztildeF-ztilde)/ztildeF
    totlossw_2 = (ztilde/ztildeF)*(muwtot+tauwtot)/Htot
    totlossw   = totlossw_1 + totlossw_2

    ! **************************************************************************
    ! Moments

    MOMPRINTNAME(1)  = "% Ajusting     "
    MOMPRINTNAME(2)  = "Average        "
    MOMPRINTNAME(3)  = "Ave. absulute  "
    MOMPRINTNAME(4)  = "Standad Des.   "
    MOMPRINTNAME(5)  = "Skewness       "
    MOMPRINTNAME(6)  = "Kurtosis       "
    MOMPRINTNAME(7)  = "% increase     "
    MOMPRINTNAME(8)  = "% abs <=5%     "
    MOMPRINTNAME(9)  = "% abs <=2.5%   "
    MOMPRINTNAME(10) = "Output losses  "
    MOMPRINTNAME(11) = "Cost errors    "
    MOMPRINTNAME(12) = "Cost setting   "
    MOMPRINTNAME(13) = "Cost timing    "
    MOMPRINTNAME(14) = "Total losses   "

    ! **************************************************************************
    ! Moments - Prices

    j1 = floor(DBLE(0.050)/(vecpchanges(2)-vecpchanges(1)))
    j2 = floor(DBLE(0.025)/(vecpchanges(2)-vecpchanges(1)))

    MOMPRINTMODEL(1)  = Avelambda*DBLE(100.0)
    MOMPRINTMODEL(2)  = AvePCh*DBLE(100.0)
    MOMPRINTMODEL(3)  = AvePChAbs*DBLE(100.0)
    MOMPRINTMODEL(4)  = StdP*DBLE(100.0)
    MOMPRINTMODEL(5)  = SkePCh
    MOMPRINTMODEL(6)  = KurPCh
    MOMPRINTMODEL(7)  = Avelambdapos*DBLE(100)
    MOMPRINTMODEL(8)  = DBLE(100)*SUM(denspchanges(nump-j1:nump+j1))
    MOMPRINTMODEL(9)  = DBLE(100)*SUM(denspchanges(nump-j2:nump+j2))
    MOMPRINTMODEL(10) = DBLE(100)*totlossp
    MOMPRINTMODEL(11) = DBLE(100)*(LossinRev_p - (atilde)/(atildeF)*((tautot+mutot)/Ntot))
    MOMPRINTMODEL(12) = DBLE(100)*(atilde)/(atildeF)*(mutot/Ntot)
    MOMPRINTMODEL(13) = DBLE(100)*(atilde)/(atildeF)*(tautot/Ntot)
    MOMPRINTMODEL(14) = DBLE(100)*totloss

    ! **************************************************************************
    ! Moments - Wages

    j1 = floor(DBLE(0.050)/(vecwchanges(2)-vecwchanges(1)))
    j2 = floor(DBLE(0.025)/(vecwchanges(2)-vecwchanges(1)))

    MOMPRINTMODEL(20) = Averho*DBLE(100)
    MOMPRINTMODEL(21) = AveWCh*DBLE(100.0)
    MOMPRINTMODEL(22) = AveWChAbs*DBLE(100.0)
    MOMPRINTMODEL(23) = StdW*DBLE(100.0)
    MOMPRINTMODEL(24) = SkeWCh
    MOMPRINTMODEL(25) = KurWCh
    MOMPRINTMODEL(26) = Averhopos*DBLE(100)
    MOMPRINTMODEL(27) = DBLE(100)*SUM(denswchanges(numw-j1:numw+j1))
    MOMPRINTMODEL(28) = DBLE(100)*SUM(denswchanges(numw-j2:numw+j2))
    MOMPRINTMODEL(29) = DBLE(100)*totlossw
    MOMPRINTMODEL(30) = DBLE(100)*(LossinRev_w - (ztilde)/(ztildeF)*((muwtot+tauwtot)/Htot))
    MOMPRINTMODEL(31) = DBLE(100)*(ztilde/ztildeF)*(muwtot)/Htot
    MOMPRINTMODEL(32) = DBLE(100)*(ztilde/ztildeF)*(tauwtot)/Htot

    ! **************************************************************************

    RETURN
  END SUBROUTINE CALCSTATS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE STATISTICS(AVE,AVEABS,STD,SKE,KUR,VAR,DIS)
    IMPLICIT NONE
    REAL(rp) , INTENT(OUT) :: AVE,AVEABS,STD,SKE,KUR
    REAL(rp) , INTENT(IN)  :: VAR(:),DIS(:)
    REAL(rp)               :: VRC
    INTEGER                :: k
    AVE = zero ; AVEABS = zero ; VRC = zero ; SKE = zero ; KUR = zero
    DO k=1,SIZE(DIS)
      AVE    = AVE    + DIS(k)*VAR(k)
      AVEABS = AVEABS + DIS(k)*ABS(VAR(k))
    END DO
    DO k=1,SIZE(DIS)
      VRC = VRC + DIS(k)*((VAR(k)-AVE)**DBLE(2.0))
      SKE = SKE + DIS(k)*((VAR(k)-AVE)**DBLE(3.0))
      KUR = KUR + DIS(k)*((VAR(k)-AVE)**DBLE(4.0))
    END DO
    STD = SQRT(VRC)
    SKE = SKE/(STD**DBLE(3.0))
    KUR = KUR/(STD**DBLE(4.0))
    RETURN
  END SUBROUTINE STATISTICS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END MODULE parameters

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
