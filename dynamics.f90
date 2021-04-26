
MODULE dynamics
IMPLICIT NONE
CONTAINS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This function takes a set of time-t distribution and value functions and a set of time-t+1
! set of fdistributions and value funtions, and computes the residuals

  FUNCTION DYNSYS(XVEC) RESULT(RESID)

    USE solution   , ONLY : ITERATE_PDIST,ITERATE_WDIST,ITERATE_L,ITERATE_V
    USE parameters , ONLY : rp,zero,one,tol,nump,numw,nums,numz,RMATRIX,RMATRIXw,&
                            delta,gamma,nu,epsilon,epsilonN,kappa_pi,kappa_lambda,lbar,S_s,S_z,&
                            p_grid,w_grid,s_grid,z_grid,mu,mu0,BIRTHDISTRIBUTION

    IMPLICIT NONE

    ! Size parameters for the dynamic system
    INTEGER , PARAMETER :: NUPV = nump*nums
    INTEGER , PARAMETER :: NUWL = numw*numz
    INTEGER , PARAMETER :: NUWZ = 5
    INTEGER , PARAMETER :: NUMB = 3
    INTEGER , PARAMETER :: NUMA = 2*NUPV + 2*NUWL + NUWZ
    INTEGER , PARAMETER :: NUMV = NUMA + NUMB
    ! Input vector with distributions and value functiosn at time t and t+1
    REAL(rp) , INTENT(IN)  :: XVEC(NUMV*2)
    ! Output vector with residuals
    REAL(rp) :: RESID(NUMA)
    ! Distributions and value functions at time t+1
    REAL(rp) :: Pdist_1(nump,nums), Wdist_1(numw,numz), V_1(nump,nums), L_1(numw,numz)
    ! Distributions and value functions at time t
    REAL(rp) :: Pdist_0(nump,nums), Wdist_0(numw,numz), V_0(nump,nums), L_0(numw,numz)
    ! Distributions and value functions at time t+1 from iterating initial distributions and value functions one period
    REAL(rp) :: Pdist_s(nump,nums), Wdist_s(numw,numz), V_s(nump,nums), L_s(numw,numz)
    ! State variables at time t
    REAL(rp) :: mlagnow,  wnow,  Infnow,  Cnow,  Nnow,  zRnow,  zSnow,  zZnow
    ! State variables at time t+1
    REAL(rp) :: mlagnext, wnext, Infnext, Cnext, Nnext, zRnext, zSnext, zZnext
    ! Other auxiliary variables
    INTEGER  :: j,iw,iz,ip,is,ips
    REAL(rp) :: T_next(nump,nump),Tw_next(numw,numw),lambda_now(nump,nums),Pi_now(nump,nums)
    REAL(rp) :: Pi_wnow(numw,numz,numw),rho_now(numw,numz),h_now(numw,numz),xp_now(numw,numz)
    REAL(rp) :: KL_pi_now,Klambda_now,Kpi_now,Delta_now,aux1
    REAL(rp) :: BirthDist0(numw,numz),mnow

    ! EXTRACTRING INFORMATION FROM INPUT VECTOR *******************************

    ! Extract value distributions
    j = 0
    DO ip=1,nump ; DO is=1,nums ; j = j + 1
      Pdist_0(ip,is) = XVEC(j) ; Pdist_1(ip,is) = XVEC(j+NUMA)
    END DO ; END DO
    DO iw=1,numw ; DO iz=1,numz ; j = j + 1
      Wdist_0(iw,iz) = XVEC(j) ; Wdist_1(iw,iz) = XVEC(j+NUMA)
    END DO ; END DO

    ! Extract state variables
    j = j + 1 ; mlagnow = XVEC(j) ; mlagnext = XVEC(j+NUMA) ; mnow = mlagnext
    j = j + 1 ; wnow    = XVEC(j) ; wnext    = XVEC(j+NUMA)
    j = j + 1 ; Infnow  = XVEC(j) ; Infnext  = XVEC(j+NUMA)

    ! Extract value functions
    DO ip=1,nump ; DO is=1,nums ; j = j + 1
      V_0(ip,is) = XVEC(j)
      V_1(ip,is) = XVEC(j+NUMA)
    END DO ; END DO
    DO iw=1,numw ; DO iz=1,numz ; j = j + 1
      L_0(iw,iz) = XVEC(j) ; L_1(iw,iz) = XVEC(j+NUMA)
    END DO ; END DO

    ! Extract jump variables
    j = j + 1 ; Cnow = XVEC(j) ; Cnext = XVEC(j+NUMA)
    j = j + 1 ; Nnow = XVEC(j) ; Nnext = XVEC(j+NUMA)

    ! Extract shocks
    zRnow = XVEC(NUMA*2+1) ; zRnext = XVEC(NUMA*2+NUMB+1)
    zSnow = XVEC(NUMA*2+2) ; zSnext = XVEC(NUMA*2+NUMB+2)
    zZnow = XVEC(NUMA*2+3) ; zZnext = XVEC(NUMA*2+NUMB+3)

    ! FIRMS PROBLEM ************************************************************
    ! Using tomorrow's value function, compute today's
    T_next = RMATRIX(Infnext)
    CALL ITERATE_V(V_s,V_1,lambda_now,Pi_now,T_next,S_s,(Cnext/Cnow)**(-gamma),Cnow,wnow)

    ! FIRMS DISTRIBUTION *******************************************************
    ! Using today's distribution, compute tomorrows's
    CALL ITERATE_PDIST(Pdist_s,Pdist_0,lambda_now,Pi_now,T_next,S_s)

    ! WORKERS PROBLEM **********************************************************
    ! Using tomorrow's value function, compute today's
    Tw_next = RMATRIXw(Infnext)
    CALL ITERATE_L(L_s,L_1,rho_now,Pi_wnow,h_now,xp_now,Tw_next,S_z,(Cnext/Cnow)**(-gamma),Cnow,wnow,Nnow)

    ! WORKERS DISTRIBUTION *****************************************************
    ! Using today's distribution, compute tomorrows's
    BirthDist0 = BIRTHDISTRIBUTION(Cnow,wnow,Nnow)
    CALL ITERATE_WDIST(Wdist_s,Wdist_0,rho_now,Pi_wnow,Tw_next,S_z,BirthDist0)

    ! LABOR IN PRICE ADJUSTMENT ************************************************
    Kpi_now = zero ; Delta_now = zero ; Klambda_now = zero ; KL_pi_now = zero
    DO ip=1,nump ; DO is = 1,nums
      KL_pi_now  = zero
      DO ips = 1,nump
        KL_pi_now = KL_pi_now + Pi_now(ips,is)*log(max(Pi_now(ips,is)*DBLE(nump),tol))
      END DO
      Kpi_now     = Kpi_now     + Pdist_0(ip,is)*lambda_now(ip,is)*KL_pi_now
      Klambda_now = Klambda_now + Pdist_0(ip,is)*(lambda_now(ip,is)*log(max(tol,lambda_now(ip,is)/lbar)) + &
                                  (one-lambda_now(ip,is))*log(max(tol,(one-lambda_now(ip,is))/(one-lbar))))
      Delta_now   = Delta_now   + Pdist_0(ip,is)*exp(s_grid(is)-epsilon*p_grid(ip))
    END DO ; END DO

    ! RESIDUALS ****************************************************************

    ! Firms distribution
    j = 0 ; DO is = 1,nums ; DO ip = 1,nump ; j = j + 1
       RESID(j) = - Pdist_1(ip,is) + Pdist_s(ip,is)
    END DO ; END DO

    ! Workers distribution
    DO iz = 1,numz ; DO iw = 1,numw ; j = j + 1
       RESID(j) = - Wdist_1(iw,iz) + Wdist_s(iw,iz)
    END DO ; END DO

    ! Money growth
    j = j + 1 ; RESID(j) = Infnow*mlagnext/mlagnow - (mu+mu0)*exp(zRnow)

    ! RFFirms' value function
    DO is = 1,nums ; DO ip = 1,nump ; j = j + 1
      RESID(j) = - V_0(ip,is) + V_s(ip,is)
    END DO ; END DO

    ! Workers' value function
    DO iz = 1,numz ; DO iw = 1,numw ; j = j + 1
      RESID(j) = - L_0(iw,iz) + L_s(iw,iz)
    END DO ; END DO

    ! Euler equation
    j = j + 1 ; RESID(j) =  one - ((one-delta)/(Cnow**(-gamma) - nu/mlagnext))*(Cnext**(-gamma)/Infnext)

    ! Aggregate labor
    j = j + 1 ; RESID(j) =  - Nnow + Delta_now*Cnow + kappa_pi*Kpi_now + kappa_lambda*Klambda_now

    ! Aggregate wage
    aux1 = zero ; DO iw = 1,numw ; DO iz = 1,numz
      aux1 = aux1 + Wdist_1(iw,iz)*exp((one-epsilonN)*(w_grid(iw)-z_grid(iz)))
    END DO ; END DO ; j = j + 1 ; RESID(j) = wnext - aux1**(one/(one-epsilonN))

    ! Aggregate price
    aux1 = zero ; DO ip = 1,nump ; DO is = 1,nums
      aux1 = aux1 + Pdist_1(ip,is)*exp(p_grid(ip)*(one-epsilon))
    END DO ; END DO ; j = j + 1 ; RESID(j) = one - aux1**(one/(one-epsilon))

    ! Checking consistency
    IF (j.NE.SIZE(RESID)) THEN
      PRINT * , '  ERROR IN DYNSYS!! RESID of inconsistent size' , j , SIZE(RESID)
    END IF

    RETURN
  END FUNCTION DYNSYS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END MODULE dynamics
