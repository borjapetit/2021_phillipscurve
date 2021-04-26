
MODULE solution
IMPLICIT NONE
CONTAINS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS SUBROUTINE COMPUTES THE VALUE FUNCTION AND DISTRIBUTION OF FIRMS

  SUBROUTINE SOLVEFIRMS
    USE parameters , ONLY : rp,tol,one,nump,nums,wbar,cbar,V,lambda,Pi,Pdist,T,S_s
    IMPLICIT NONE
    INTEGER  :: iter
    REAL(rp) :: test,V0(nump,nums),V1(nump,nums),Dist1(nump,nums),Dist0(nump,nums)

    ! ----------------------------------------------------------------------------------
    ! Value function:
    iter = 0 ; test = 99.999999 ; V1 = V ; V0 = V
    DO WHILE (ABS(test).gt.tol .and. iter.lt.15000) ; iter = iter + 1
      CALL ITERATE_V(V0,V1,lambda,Pi,T,S_s,one,cbar,wbar)
      test = MAXVAL(ABS(V1-V0))
      V1 = V0
    END DO
    V = V1

    ! ----------------------------------------------------------------------------------
    ! Distribution:
    iter = 0 ; test = 99.999 ; Dist0 = Pdist ; Dist1 = Pdist
    DO WHILE (abs(test).gt.tol .and. iter.lt.15000) ; iter = iter + 1
      CALL ITERATE_PDIST(Dist1,Dist0,lambda,Pi,T,S_s)
      test = MAXVAL(ABS(Dist1-Dist0))
      Dist0 = Dist1/SUM(Dist1)
    END DO
    Pdist = Dist1/SUM(Dist1)

    RETURN
  END SUBROUTINE SOLVEFIRMS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS SUBROUTINE COMPUTES THE VALUE FUNCTION AND DISTRIBUTION OF WORKERS

  SUBROUTINE SOLVEWORKERS
    USE parameters , ONLY : BIRTHDISTRIBUTION,rp,tol,one,numw,numz,wbar,nbar,cbar,L,H,rho,Pi_w,Wdist,Xp,Tw,S_z
    IMPLICIT NONE
    INTEGER  :: iter
    REAL(rp) :: test,L0(numw,numz),L1(numw,numz),Dist0(numw,numz),Dist1(numw,numz),Bdist(numw,numz)

    ! ----------------------------------------------------------------------------------
    ! Value function:
    iter = 0 ; test = 99.999 ; L1 = L ; L0 = L
    DO WHILE (test.gt.tol .and. iter.lt.15000) ; iter = iter + 1
      CALL ITERATE_L(L0,L1,rho,Pi_w,H,Xp,Tw,S_z,one,cbar,wbar,nbar)
      test = MAXVAL(ABS(L1-L0))
      L1 = L0
    END DO
    L = L1

    ! ----------------------------------------------------------------------------------
    ! Distribution:
    Bdist = BIRTHDISTRIBUTION(cbar,wbar,nbar)
    iter  = 0 ; test = 99.999 ; Dist0 = Wdist ; Dist1 = Wdist
    DO WHILE (test.gt.tol .and. iter.lt.15000) ; iter = iter + 1
      CALL ITERATE_WDIST(Dist1,Dist0,rho,Pi_w,Tw,S_z,Bdist)
      test = MAXVAL(ABS(Dist1-Dist0))
      Dist0 = Dist1/SUM(Dist1)
    END DO
    Wdist = Dist1/SUM(Dist1)

    RETURN
  END SUBROUTINE SOLVEWORKERS

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS SUBROUTINE TAKES TOMORROWS VALUE FUNCTION OF FIRMS, AND COMPUTES TODAY'S ONE

  SUBROUTINE ITERATE_V(VTODAY,VTOMORROW,LAMBDAMAT,PIMAT,TMAT,SMAT,DISC,CONS,WAGE)

    USE toolkit    , ONLY : INTERPOLATION
    USE parameters , ONLY : rp,tol,zero,one,nump,nums,lbar,kappa_lambda,&
                            kappa_pi,p_grid,s_grid,delta,epsilon,LABUS,Lmu
    IMPLICIT NONE

    REAL(rp) , INTENT(IN)  :: VTOMORROW(nump,nums),TMAT(nump,nump),SMAT(nums,nums)
    REAL(rp) , INTENT(IN)  :: DISC,CONS,WAGE
    REAL(rp) , INTENT(OUT) :: VTODAY(nump,nums),LAMBDAMAT(nump,nums),PIMAT(nump,nums)

    INTEGER  :: ip,is
    REAL(rp) :: mat(nump,nums),vecp(nump),Ve(nump,nums),EV,Vtilde,KL_pi,KL_lambda

    Ve        = zero
    PIMAT     = zero
    LAMBDAMAT = zero
    VTODAY    = zero
    DO is = 1,nums ; DO ip = 1,nump
      mat(ip,is) = SUM(TMAT(:,ip)*VTOMORROW(:,is))
    END DO ; END DO
    DO is = 1,nums ; DO ip = 1,nump
      Ve(ip,is) = (one-delta)*DISC*SUM(mat(ip,:)*SMAT(:,is))
    END DO ; END DO
    DO is = 1,nums
      DO ip = 1,nump
        vecp(ip) = (one/DBLE(nump))*exp((Ve(ip,is)-maxval(Ve(:,is)))/(kappa_pi*WAGE))
      END DO
      DO ip = 1,nump
        PIMAT(ip,is) = vecp(ip)/sum(vecp)
      END DO
      EV = zero ; KL_pi = zero
      DO ip = 1,nump
        EV    = EV    + PIMAT(ip,is)*Ve(ip,is)
        KL_pi = KL_pi + PIMAT(ip,is)*LOG(MAX(PIMAT(ip,is)*DBLE(nump),tol))
      END DO
      Vtilde = EV - WAGE*kappa_pi*KL_pi
      DO ip = 1,nump
        LAMBDAMAT(ip,is) = lbar/( lbar + (one-lbar)*exp(-(Vtilde-Ve(ip,is))/(WAGE*kappa_lambda)) )
        KL_lambda        = LAMBDAMAT(ip,is)*log(max(tol,LAMBDAMAT(ip,is)/lbar)) + &
                           (one-LAMBDAMAT(ip,is))*log(max(tol,(one-LAMBDAMAT(ip,is))/(one-lbar)) )
        LABUS(ip,is)     = LAMBDAMAT(ip,is)*kappa_pi*KL_pi + kappa_lambda*KL_lambda
        LMU(ip,is)       = kappa_lambda*KL_lambda
        VTODAY(ip,is)    = ( CONS*exp(p_grid(ip)*(one-epsilon)) - WAGE*CONS*exp(-p_grid(ip)*epsilon)*exp(s_grid(is)) ) + &
                           Ve(ip,is) + LAMBDAMAT(ip,is)*( EV - WAGE*kappa_pi*KL_pi - Ve(ip,is) ) - WAGE*kappa_lambda*KL_lambda
      END DO
    END DO

    RETURN
  END SUBROUTINE ITERATE_V

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS SUBROUTINE TAKES TODAY'S DISTRIBUTION OF FIRMS, AND COMPUTES TOMORROWS'S ONE

  SUBROUTINE ITERATE_PDIST(PNEW,PDOLD,LAMBDAMAT,PIMAT,TMAT,SMAT)
    USE toolkit    , ONLY : INTERPOLATION
    USE parameters , ONLY : rp,zero,one,nump,nums
    IMPLICIT NONE
    REAL(rp) , INTENT(IN)  :: PDOLD(nump,nums),LAMBDAMAT(nump,nums),PIMAT(nump,nums)
    REAL(rp) , INTENT(IN)  :: TMAT(nump,nump),SMAT(nums,nums)
    REAL(rp) , INTENT(OUT) :: PNEW(nump,nums)
    REAL(rp)               :: Dist1(nump,nums),Dist2(nump,nums),Dist3(nump,nums)
    INTEGER                :: ip,is
    Dist1 = zero
    DO is = 1,nums ; DO ip = 1,nump
      Dist1(ip,is) = (one-LAMBDAMAT(ip,is))*PDOLD(ip,is) + PIMAT(ip,is)*SUM(LAMBDAMAT(:,is)*PDOLD(:,is))
    END DO ; END DO
    DO is = 1,nums ; DO ip = 1,nump
      Dist2(ip,is) = SUM(TMAT(ip,:)*Dist1(:,is))
    END DO ; END DO
    DO is = 1,nums ; DO ip = 1,nump
      Dist3(ip,is) = SUM(Dist2(ip,:)*SMAT(is,:))
    END DO ; END DO
    PNEW = Dist3
    RETURN
  END SUBROUTINE ITERATE_PDIST

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS SUBROUTINE TAKES TOMORROWS VALUE FUNCTION OF WORKERS, AND COMPUTES TODAY'S ONE

  SUBROUTINE ITERATE_L(LTODAY,LTOMORROW,RHOMAT,PIMAT,HMAT,XMAT,TMAT,SMAT,DISC,CONS,WAGE,NTOT)

    USE toolkit    , ONLY : INTERPOLATION
    USE parameters , ONLY : rp,numthreads,tol,zero,one,numw,numz,rhobar,kappa_w,Ld,Lrho,&
                            kappa_rho,w_grid,z_grid,delta,gamma,labelast,chi,epsilonN
    IMPLICIT NONE

    REAL(rp) , INTENT(IN)  :: LTOMORROW(numw,numz),TMAT(numw,numw),SMAT(numz,numz)
    REAL(rp) , INTENT(IN)  :: DISC,CONS,WAGE,NTOT
    REAL(rp) , INTENT(OUT) :: LTODAY(numw,numz),RHOMAT(numw,numz),PIMAT(numw,numz,numw),XMAT(numw,numz),HMAT(numw,numz)

    INTEGER  :: iw,iz,iws,iterbeta
    REAL(rp) :: Le(numw,numz),mat(numw,numz)
    REAL(rp) :: testbeta,xpnow,rhowz,KL_w,piwz(numw),labnow
    REAL(rp) :: expG,KL_rho,vecw(numw),maxvecw,maxlab,minlab

    DO iw = 1,numw ; DO iz = 1,numz ; Ld(iw,iz) = zero ; Lrho(iw,iz) = zero
      Ld(iw,iz) = (WAGE**epsilonN)*NTOT*exp( z_grid(iz)*(epsilonN-one) - w_grid(iw)*epsilonN )
    END DO ; END DO
    DO iw = 1,numw ; DO iz = 1,numz
      mat(iw,iz) = SUM(TMAT(:,iw)*LTOMORROW(:,iz))
    END DO ; END DO
    DO iw = 1,numw ; DO iz = 1,numz
      Le(iw,iz) = (one-delta)*DISC*SUM(mat(iw,:)*SMAT(:,iz))
    END DO ; END DO
    !$OMP PARALLEL NUM_THREADS(numthreads), DEFAULT(SHARED), &
    !$OMP PRIVATE(labnow,maxlab,minlab,xpnow,iws,vecw,maxvecw,piwz,KL_w,expG,rhowz,KL_rho,testbeta,iterbeta)
    !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
    DO iw = 1,numw ; DO iz = 1,numz
      maxlab   = 10.00000
      minlab   = 0.000000
      iterbeta = 0
      testbeta = 99.9999
      DO WHILE ( abs(testbeta).gt.tol .and. iterbeta.lt.10000  .and. abs(maxlab-minlab).ge.tol ) ; iterbeta = iterbeta + 1
        labnow = DBLE(0.5)*( maxlab + minlab )
        xpnow  = ( chi*(labnow**labelast) )*(CONS**gamma)
        DO iws = 1,numw
          vecw(iws) = Le(iws,iz)/(xpnow*kappa_w)
        END DO ; maxvecw = maxval(vecw)
        DO iws = 1,numw
          vecw(iws) = exp(vecw(iws) - maxvecw)
        END DO ; maxvecw = sum(vecw)
        DO iws = 1,numw
          piwz(iws) = vecw(iws)/maxvecw
        END DO
        KL_w  = zero
        DO iws = 1,numw
          KL_w = KL_w + piwz(iws)*LOG(MAX(tol*DBLE(numw),piwz(iws)*DBLE(numw)))
        END DO
        expG     = SUM(piwz(:)*Le(:,iz)) - Le(iw,iz) - kappa_w*xpnow*KL_w
        rhowz    = rhobar/(rhobar + (one-rhobar)*EXP(-expG/(kappa_rho*xpnow)))
        KL_rho   = rhowz*LOG(MAX(tol,rhowz/rhobar)) + (one-rhowz)*LOG(MAX(tol,(one-rhowz)/(one-rhobar)))
        testbeta = labnow - Ld(iw,iz) - kappa_w*rhowz*KL_w - kappa_rho*KL_rho
        IF ( testbeta.gt.zero ) maxlab = labnow
        IF ( testbeta.lt.zero ) minlab = labnow
      END DO
      HMAT(iw,iz)    = Ld(iw,iz) + kappa_w*rhowz*KL_w + kappa_rho*KL_rho
      XMAT(iw,iz)    = xpnow
      PIMAT(iw,iz,:) = piwz
      RHOMAT(iw,iz)  = rhowz
      LRHO(iw,iz)    = kappa_rho*KL_rho
      LTODAY(iw,iz)  = EXP(w_grid(iw))*Ld(iw,iz) - &
                         chi*( (HMAT(iw,iz)**(one+labelast))/(one+labelast) )*(CONS**gamma) + &
                         (one-RHOMAT(iw,iz))*Le(iw,iz) + RHOMAT(iw,iz)*SUM(PIMAT(iw,iz,:)*Le(:,iz))
    END DO ; END DO
    !$OMP END DO
    !$OMP END PARALLEL

    RETURN
  END SUBROUTINE ITERATE_L

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS SUBROUTINE TAKES TODAY'S DISTRIBUTION OF WORKERS, AND COMPUTES TOMORROWS'S ONE

  SUBROUTINE ITERATE_WDIST(WNEW,WDOLD,RHOMAT,PIMAT,TMAT,SMAT,BIRTHMAT)
    USE toolkit    , ONLY : INTERPOLATION
    USE parameters , ONLY : rp,zero,one,numw,numz,deathprob
    IMPLICIT NONE
    REAL(rp) , INTENT(IN)  :: WDOLD(numw,numz),RHOMAT(numw,numz),PIMAT(numw,numz,numw)
    REAL(rp) , INTENT(IN)  :: TMAT(numw,numw),SMAT(numz,numz),BIRTHMAT(numw,numz)
    REAL(rp) , INTENT(OUT) :: WNEW(numw,numz)
    REAL(rp)               :: Dist1(numw,numz),Dist2(numw,numz),Dist3(numw,numz)
    INTEGER                :: iw,iz,iws
    Dist1 = zero
    Dist2 = zero
    DO iz  = 1,numz ; DO iw = 1,numw
      Dist1(iw,iz) = (one-RHOMAT(iw,iz))*WDOLD(iw,iz)
    END DO ; END DO
    DO iws = 1,numw ; DO iz  = 1,numz ; DO iw = 1,numw
      Dist2(iws,iz) = Dist2(iws,iz) + RHOMAT(iw,iz)*WDOLD(iw,iz)*PIMAT(iw,iz,iws)
    END DO ; END DO ; END DO
    Dist3 = Dist1 + Dist2
    Dist1 = zero
    Dist2 = zero
    DO iz = 1,numz ; DO iw = 1,numw
      Dist1(iw,iz) = SUM(TMAT(iw,:)*Dist3(:,iz))
    END DO ; END DO
    DO iz = 1,numz ; DO iw = 1,numw
      Dist2(iw,iz) = SUM(Dist1(iw,:)*SMAT(iz,:))
    END DO ; END DO
    Dist3 = zero
    DO iz = 1,numz ; DO iw = 1,numw
      Dist3(iw,iz) = (one-deathprob)*Dist2(iw,iz) + deathprob*BIRTHMAT(iw,iz)
    END DO ; END DO
    WNEW = Dist3
    RETURN
  END SUBROUTINE ITERATE_WDIST

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

END MODULE solution
