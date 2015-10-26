!############################################################################
SUBROUTINE VINIT(Nset, eam, param, oxygn, sys, q, cell, psi, dt)
USE DATAFMT
!
TYPE(AMNT):: Nset
TYPE(METL):: eam(Nmetal)
TYPE(POTN):: param(Natom)
TYPE(OXID):: oxygn
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
TYPE(SORT):: cell(Nset%Ncell_all)
REAL(DP)  :: psi(Nset%fft(1),Nset%fft(2),Nset%fft(3)), dt
!
TYPE(PTCL):: old(Nset%Npt_all)
TYPE(CORR):: pair
INTEGER(GI):: i, j, id, Npt(Natom)
REAL(DP):: xm, lambda(Natom), T0, xv(Nset%Npt_all,3), mv2(Natom), sumv(3), rand
!
! Velocity initialization by random number - rigid motion is removed
sumv(:) = 0.0D0
DO i = 1, Nset%Npt_all
   DO j=1, 3
      CALL RANDOM_NUMBER(rand)
      xv(i,j) = rand - 0.5D0
      sumv(j) = sumv(j) + xv(i,j)
   END DO
END DO
!
sumv(:) = sumv(:)/DBLE(Nset%Npt_all)
mv2(:) = 0.0D0
Npt(:) = 0
DO i=1, Nset%Npt_all
   xv(i,:)= xv(i,:) - sumv(:)
   id = q(i)%id
   xm = param(id)%xm   
   mv2(id) = mv2(id) + xm*(xv(i,1)**2 + xv(i,2)**2 + xv(i,3)**2)
   Npt(id) = Npt(id) + 1
END DO
!
lambda(:) = 0.0D0
DO i=1, Natom
   IF (Npt(i) > 0) THEN
      T0 = mv2(i)/DBLE(Npt(i)*3)
      lambda(i) = DSQRT(sys%T_given/T0)
   END IF
END DO
!
DO i = 1, Nset%Npt_all
   id = q(i)%id
   q(i)%xv(:) = xv(i,:)*lambda(id)
   q(i)%ff(:) = 0.0D0
   old(i)%xx(:) = q(i)%xx(:)
   old(i)%xv(:) = q(i)%xv(:)
END DO
!
! Initial force calculation - compensate the motion by the initial velocity
CALL VV_NVE1(Nset, param, q, -dt)
pair%dr = sys%Rcut/DBLE(Nrdf)
pair%rdf(:,:,:) = 0
CALL FORCE(Nset, eam, oxygn, sys, param, pair, q, cell)
CALL PME(Nset, sys, q, psi)
sys%mv2 = 0.0D0
DO i=1, Nset%Npt_all
   id = q(i)%id
   q(i)%xx(:) = old(i)%xx(:)
   q(i)%xv(:) = old(i)%xv(:)
   sys%mv2 = sys%mv2 + xm*(q(i)%xv(1)**2 + q(i)%xv(2)**2 + q(i)%xv(3)**2 )
END DO
!
RETURN
END SUBROUTINE Vinit
!
!############################################################################
SUBROUTINE VV_NVE1(Nset, param, q, dt)
USE DATAFMT
!
TYPE(AMNT):: Nset
TYPE(POTN):: param(Natom)
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, id
REAL(DP):: xm, x(3), box(3), box_z
!
box(:) = Nset%box(:)
box_z = box(3)/2.D0
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%xm
   q(i)%xv(:) = q(i)%xv(:) + 0.5D0 * dt * q(i)%ff(:)/xm
   x(:) = q(i)%xx(:) + dt*q(i)%xv(:)
   q(i)%xx(1:3) = x(1:3) - box(1:3)*DNINT(x(1:3)/box(1:3))
END DO
!
RETURN
END SUBROUTINE VV_NVE1
!
!
!############################################################################
SUBROUTINE VV_NVE1_REFLECT(Nset, param, q, dt)
USE DATAFMT
!
TYPE(AMNT):: Nset
TYPE(POTN):: param(Natom)
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, id
REAL(DP):: xm, x(3), box(3), box_z
!
box(:) = Nset%box(:)
box_z = box(3)/2.D0
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%xm
   q(i)%xv(:) = q(i)%xv(:) + 0.5D0 * dt * q(i)%ff(:)/xm
   x(:) = q(i)%xx(:) + dt*q(i)%xv(:)
   q(i)%xx(1:2) = x(1:2) - box(1:2)*DNINT(x(1:2)/box(1:2))
   !
   ! Reflecting boundary condition on the top
   IF (x(3) > box_z) THEN
      x(3) = box(3) - x(3)
      q(i)%xv(3) = -DABS(q(i)%xv(3))
   END IF
   q(i)%xx(3) = x(3)
END DO
!
RETURN
END SUBROUTINE VV_NVE1_REFLECT
!
!############################################################################
SUBROUTINE VV_NVE2(Nset, param, sys, q, dt)
USE DATAFMT
!
TYPE(AMNT):: Nset
TYPE(POTN):: param(Natom)
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, k, id
REAL(DP)   :: xm, mv2_tmp
!
sys%mv2 = 0.0D0
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%xm
   q(i)%xv(:) = q(i)%xv(:) + 0.5D0 * dt * q(i)%ff(:)/xm
   mv2_tmp = 0.0D0
   DO k=1,3
      mv2_tmp = mv2_tmp + q(i)%xv(k)**2
   END DO
   sys%mv2 = sys%mv2 + xm*mv2_tmp
END DO
!
RETURN
END SUBROUTINE VV_NVE2
!
!############################################################################
SUBROUTINE VV_NVT2_BEREND(Nset, param, sys, q, dt)
USE DATAFMT
!
TYPE(AMNT):: Nset
TYPE(POTN):: param(Natom)
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, k, id, Npt(Natom)
REAL(DP):: xm, mv2_tmp, T0, lambda(Natom), mv2_sum(Natom)
REAL(DP), PARAMETER:: C_BEREND = 0.025D0
!
!
mv2_sum(:) = 0.0D0; lambda(:) = 0.0D0
Npt(:) = 0
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%xm
   q(i)%xv(:) = q(i)%xv(:) + 0.5D0 * dt * q(i)%ff(:)/xm
   mv2_tmp = 0.0D0
   DO k=1,3
      mv2_tmp = mv2_tmp + xm*q(i)%xv(k)**2
   END DO
   mv2_sum(id) = mv2_sum(id) + mv2_tmp
   Npt(id) = Npt(id) + 1
END DO
!
DO i=1, Natom
   IF (Npt(i)> 0) THEN
      T0 = mv2_sum(i)/DBLE(3*Npt(i))
      lambda(i) = DSQRT(1.D0 + C_BEREND*(sys%T_given/T0 - 1.0D0))
   END IF
END DO
!
sys%mv2 = 0.0D0
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%xm
   q(i)%xv(:) = q(i)%xv(:)*lambda(id)
   mv2_tmp = 0.0D0
   DO k=1,3
      mv2_tmp = mv2_tmp + xm*q(i)%xv(k)**2
   END DO
   sys%mv2 = sys%mv2 + mv2_tmp
END DO
!
RETURN
END SUBROUTINE VV_NVT2_BEREND
!
!############################################################################
SUBROUTINE REMOVE_RIGID_MOTION(Nset, q)
USE DATAFMT
!
TYPE(AMNT):: Nset
TYPE(PTCL):: q(Nset%Npt_all)
!
INTEGER(GI):: i
REAL   (DP):: sumv(3), avg_z
!
sumv(:) = 0.0D0
avg_z = 0.0D0
DO i=1, Nset%Npt_all
   sumv(:) = sumv(:) + q(i)%xv(:)
   avg_z   = avg_z   + q(i)%xx(3)
END DO
sumv(:) = sumv(:)/DBLE(Nset%Npt_all)
avg_z   = avg_z  /DBLE(Nset%Npt_all)
DO i=1, Nset%Npt_all
   q(i)%xv(:) = q(i)%xv(:) - sumv(:)
   q(i)%xx(3) = q(i)%xx(3) - avg_z
END DO
!
RETURN
END SUBROUTINE REMOVE_RIGID_MOTION
!
!############################################################################
SUBROUTINE VV_NVT1_STOCH(Nset, param, q, dt)
USE DATAFMT
!
TYPE(AMNT):: Nset
TYPE(POTN):: param(Natom)
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, id
REAL(DP):: xm, x(3), box(3), alpha
!
alpha = Nset%alpha
box(:) = Nset%box(:)
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%xm
   q(i)%xv(:) = q(i)%xv(:)*(1.D0 - 0.5D0*alpha*dt/xm) + &
        & 0.5D0 * dt * q(i)%ff(:)/xm
   x(:) = q(i)%xx(:) + dt*q(i)%xv(:)
   q(i)%xx(:) = x(:) - box(:)*DNINT(x(:)/box(:))
END DO
!
RETURN
END SUBROUTINE VV_NVT1_STOCH
!
!############################################################################
SUBROUTINE VV_NVT2_STOCH(Nset, param, sys, q, dt)
USE DATAFMT
!
TYPE(AMNT):: Nset
TYPE(POTN):: param(Natom)
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: dt
!
INTEGER(GI):: i, k, id
REAL(DP):: xm, mv2_tmp, fluct, sigma, eta
!
!
alpha = Nset%alpha
sigma = 1.0D0
beta = DSQRT(2.D0*alpha*sys%T_given/dt)
sys%mv2 = 0.0D0
DO i=1, Nset%Npt_all
   id = q(i)%id
   xm = param(id)%xm
   mv2_tmp = 0.0D0
   DO k=1,3
      eta = fluct(sigma)
      q(i)%ff(k) = q(i)%ff(k) + eta*beta
      q(i)%xv(k) = (q(i)%xv(k) + 0.5D0*dt*q(i)%ff(k)/xm) / &
           &       (1.D0 + 0.5D0*alpha*dt/xm)
      mv2_tmp = mv2_tmp + xm*q(i)%xv(k)**2
   END DO
   sys%mv2 = sys%mv2 + mv2_tmp
END DO
!
!
RETURN
END SUBROUTINE VV_NVT2_STOCH
