!
!##############################################################################
SUBROUTINE CTIP_SOLVER(Nset, param, sys, q, cell, psi)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(POTN):: param(Natom)
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
TYPE(SORT):: cell(Nset%Ncell_all)
REAL(DP)  :: psi(Nset%fft(1),Nset%fft(2),Nset%fft(3))

!
REAL(DP):: XJZ(Nset%Npt_all,3), q_i(Nset%Npt_all), g_sum, g_sum_old, d_sum
REAL(DP):: d(Nset%Npt_all), p(Nset%Npt_all), g(Nset%Npt_all), Xall, q_sum
REAL(DP):: afb_i, xi_i, EXP_rxi_i, EXP_rxi_i_r, a, ar, box(3)
REAL(DP):: afb_j, xi_j, EXP_rxi_j, EXP_rxi_j_r, dx(3), rx(3), r2, r, Esum, dq
REAL(DP):: xir, fifj, kc_const, alpha, V, Eself, c_pme, tmp, rcut2
REAL(DP):: asym_mat(Nset%Npt_all, Nnb_max), phi(Nset%Npt_all)
REAL(DP):: sym_mat(Nset%Npt_all,Nnb_max), qmx(Nset%Npt_all, 2)
REAL(DP):: u(3), z(3), AX(4,3), Qlocal(Nset%Npt_all,4,4,4), xfft(3), Error
INTEGER(GI):: Npt_all, Ncell_all, id_i, id_j, grid(Nset%Npt_all,4,4,4,3)
INTEGER(GI):: fft(3), nn(3), list(Nset%Npt_all,Nnb_max), nlist(Nset%Npt_all)
INTEGER(GI):: i, j, k, l, m, n, id_pair, Npt_pair, Npt_local, l_i, l_j
REAL(DP), PARAMETER:: E_crit = 1.D-3
!
! Preparation of PME data array
Npt_all = Nset%Npt_all; Xall = DBLE(Npt_all)
Ncell_all = Nset%Ncell_all
V = Nset%V; a = sys%a; kc_const = kc*2.D0*a/SQRTPI; rcut2 = sys%Rcut2
box(:) = Nset%box(:); fft(:) = Nset%fft(:); xfft(:) = DBLE(fft(:))
DO i=1, Npt_all
   id_i = q(i)%id
   q_i(i) = q(i)%q
   XJZ(i,1) = param(id_i)%chi
   XJZ(i,2) = param(id_i)%J - kc_const
   XJZ(i,3) = param(id_i)%Z
   qmx(i,1) = param(id_i)%qmin; qmx(i,2) = param(id_i)%qmax
   !
   u(:) = xfft(:)*(q(i)%xx(:) + box(:)*0.5D0)/box(:)
   z(:) = DINT(u(:)) - u(:) + 1.D0
   AX(1,:)  =       z(:)**3                                      /6.D0
   z(:) = z(:) + 1.D0
   AX(2,:) = (-3.D0*z(:)**3 + 12.D0*z(:)**2 - 12.D0*z(:) +  4.D0)/6.D0
   z(:) = z(:) + 1.D0
   AX(3,:) = ( 3.D0*z(:)**3 - 24.D0*z(:)**2 + 60.D0*z(:) - 44.D0)/6.D0
   z(:) = z(:) + 1.D0
   AX(4,:) = (     -z(:)**3 + 12.D0*z(:)**2 - 48.D0*z(:) + 64.D0)/6.D0
   DO j=1,4; nn(1) = j
      DO k=1,4; nn(2) = k
         DO l=1,4; nn(3) = l
            grid(i,j,k,l,:) = MOD(INT(u(:)) + nn(:) - 2 + fft(:), fft(:)) + 1
            Qlocal(i,j,k,l) = AX(j,1)*AX(k,2)*AX(l,3)
         END DO
      END DO
   END DO
END DO
c_pme = kc * xfft(1)*xfft(2)*xfft(3)/V/PI
!
! Preparation of Coulomb integration data
! asym_mat(i,j) = [j|fi] - [fi|fj]
! sym_mat(i,j) = [fi|fj] - 1/rij + erfc(arij)/rij
nlist(:) = 0
list(:,:) = 0
asym_mat(:,:) = 0.0D0; sym_mat(:,:) = 0.0D0
DO m=1, Ncell_all
   !
   ! Interactions inside of a single cell
   Npt_local = cell(m)%Npt
   DO i=1, Npt_local-1
      l_i = cell(m)%link(i)
      id_i = q(l_i)%id
      xi_i = param(id_i)%xi
      DO j=i+1, Npt_local
         l_j = cell(m)%link(j)
         dx(:) = q(l_i)%xx(:) - q(l_j)%xx(:)
         rx(:) = dx(:) - box(:)*DNINT(dx(:)/box(:))
         r2 = rx(1)**2 + rx(2)**2 + rx(3)**2
         IF (r2 < rcut2) THEN
            r = DSQRT(r2)
            ar = a*r
            id_j = q(l_j)%id
            xi_j = param(id_j)%xi
            EXP_rxi_j = DEXP(-2.D0*xi_j*r); EXP_rxi_j_r = EXP_rxi_j/r
            afb_j = -xi_j*EXP_rxi_j - EXP_rxi_j_r
            IF( id_i == id_j) THEN
               EXP_rxi_i = EXP_rxi_j; EXP_rxi_i_r = EXP_rxi_j_r
               afb_i = afb_j
               xir = xi_i*r
               fifj = - EXP_rxi_j*(1.D0 + xir*(param(id_i)%fafb(id_i,2) + &
                    & xir*(param(id_i)%fafb(id_i,3) + &
                    & xir*param(id_i)%fafb(id_i,4))))/r
            ELSE
               EXP_rxi_i = DEXP(-2.D0*xi_i*r); EXP_rxi_i_r = EXP_rxi_i/r
               afb_i = -xi_i*EXP_rxi_i - EXP_rxi_i_r
               fifj = - EXP_rxi_i* &
                    & (param(id_i)%fafb(id_j,1) + param(id_i)%fafb(id_j,3)/r) &
                    & - EXP_rxi_j* &
                    & (param(id_i)%fafb(id_j,2) + param(id_i)%fafb(id_j,4)/r)
            END IF
            nlist(l_i) = nlist(l_i) + 1
            nlist(l_j) = nlist(l_j) + 1
            list(l_i,nlist(l_i)) = l_j
            list(l_j,nlist(l_j)) = l_i
            asym_mat(l_i,nlist(l_i)) = kc*(afb_i - fifj)
            asym_mat(l_j,nlist(l_j)) = kc*(afb_j - fifj)
            tmp = kc*(fifj + DERFC(ar)/r)
            sym_mat(l_i,nlist(l_i)) = tmp
            sym_mat(l_j,nlist(l_j)) = tmp
         END IF
      END DO
   END DO
   !
   ! Interactions nearby cell (13 pairs)
   DO n=1, 13
      id_pair = cell(m)%pair(n)
      Npt_pair = cell(id_pair)%Npt
      DO j=1, Npt_pair
         l_j = cell(id_pair)%link(j)
         id_j = q(l_j)%id
         xi_j = param(id_j)%xi
         DO i=1, Npt_local
            l_i = cell(m)%link(i)
            dx(:) = q(l_i)%xx(:) - q(l_j)%xx(:)
            rx(:) = dx(:) - box(:)*DNINT(dx(:)/box(:))
            r2 = rx(1)**2 + rx(2)**2 + rx(3)**2
            IF (r2 < rcut2) THEN
               r = DSQRT(r2)
               id_i = q(l_i)%id
               xi_i = param(id_i)%xi
               ar = a*r
               EXP_rxi_j = DEXP(-2.D0*xi_j*r); EXP_rxi_j_r = EXP_rxi_j/r
               afb_j = -xi_j*EXP_rxi_j - EXP_rxi_j_r
               IF( id_i == id_j) THEN
                  EXP_rxi_i = EXP_rxi_j; EXP_rxi_i_r = EXP_rxi_j_r
                  afb_i = afb_j
                  xir = xi_i*r
                  fifj = - EXP_rxi_j*(1.D0 + xir*(param(id_i)%fafb(id_i,2) + &
                       & xir*(param(id_i)%fafb(id_i,3) + &
                       & xir*param(id_i)%fafb(id_i,4))))/r
               ELSE
                  EXP_rxi_i = DEXP(-2.D0*xi_i*r); EXP_rxi_i_r = EXP_rxi_i/r
                  afb_i = -xi_i*EXP_rxi_i - EXP_rxi_i_r
                  fifj = - EXP_rxi_i* &
                       & (param(id_i)%fafb(id_j,1)+param(id_i)%fafb(id_j,3)/r)&
                       & - EXP_rxi_j* &
                       & (param(id_i)%fafb(id_j,2)+param(id_i)%fafb(id_j,4)/r)
               END IF
               nlist(l_i) = nlist(l_i) + 1
               nlist(l_j) = nlist(l_j) + 1
               IF (nlist(l_i) > Nnb_max .OR. nlist(l_j) > Nnb_max) &
                    & STOP "=== neighbor particle limit crashed !!! === "
               list(l_i,nlist(l_i)) = l_j
               list(l_j,nlist(l_j)) = l_i
               asym_mat(l_i,nlist(l_i)) = kc*(afb_i - fifj)
               asym_mat(l_j,nlist(l_j)) = kc*(afb_j - fifj)
               tmp = kc*(fifj + DERFC(ar)/r)
               sym_mat(l_i,nlist(l_i)) = tmp
               sym_mat(l_j,nlist(l_j)) = tmp
            END IF
         END DO
      END DO
   END DO
END DO
!
! Iterative conjugate gradient solver
Error = 1.0; p(:) = 0.0D0; g_sum = 1.0D0; k = 0
DO WHILE (Error > E_crit)
   CALL PHI_MESH(Npt_all, fft, psi, Qlocal, q_i, grid, phi); phi(:)=phi(:)*c_pme
   CALL dUdq(Npt_all, XJZ, asym_mat, sym_mat, list, nlist, q_i, qmx, phi, g)
   g_sum_old = g_sum; g_sum = 0.0D0
   DO i=1, Npt_all
      g_sum = g_sum + g(i)
   END DO
   d(:) = g(:) + p(:)*g_sum/g_sum_old
   d_sum = 0.0D0
   DO i=1, Npt_all
      d_sum = d_sum + d(i)
   END DO; d_sum = d_sum/Xall
   p(:) = d(:) - d_sum
   !
   CALL BRENT_SOLVER(Npt_all, fft, psi, Qlocal, grid, c_pme, XJZ, asym_mat, &
        & sym_mat, list, nlist, q_i, qmx, p, alpha)
   q_i(:) = q_i(:) + alpha*p(:)
   !   
   Error= MAXVAL(DABS(p(:)))*DABS(alpha)
   !WRITE(*,100) Error, k
   k = k + 1
END DO
!WRITE(*,100) Error, k
!100 FORMAT("# CG error =", ES9.2, " at ", I3, "th loop")
!
Eself = 0.0D0
q_sum = 0.0D0
Esum = 0.0D0
DO i=1, Npt_all
   q(i)%q = q_i(i)
   Eself = Eself - q(i)%q**2
   q_sum = q_sum + q_i(i)
   Esum = Esum + XJZ(i,1)*q_i(i) + 0.5D0*XJZ(i,2)*q_i(i)**2
   dq = q_i(i) - qmx(i,1)
   IF (dq < 0.0D0) Esum = Esum + w*dq**2
   dq = q_i(i) - qmx(i,2)
   IF (dq > 0.0D0) Esum = Esum + w*dq**2
END DO
!WRITE(*,150) q_sum
!150 FORMAT("# Charge sum of CTIP loop = ", ES11.4)
sys%Eself = Eself*kc*a/SQRTPI + Esum
!
RETURN
END SUBROUTINE CTIP_SOLVER
!
! Brent solver
!##############################################################################
SUBROUTINE BRENT_SOLVER(Npt_all, fft, psi, Qlocal, grid, c_pme, XJZ, asym_mat, &
     & sym_mat, list, nlist, q_i, qmx, p, alpha)
USE DATAFMT
IMPLICIT NONE
INTEGER(GI):: Npt_all, fft(3)
REAL   (DP):: psi(fft(1),fft(2),fft(3)), Qlocal(Npt_all,4,4,4), c_pme
INTEGER(GI):: grid(Npt_all,4,4,4,3), list(Npt_all, Nnb_max), nlist(Npt_all)
REAL   (DP):: XJZ(Npt_all,3), asym_mat(Npt_all,Nnb_max), p(Npt_all)
REAL   (DP):: q_i(Npt_all), qmx(Npt_all,2), sym_mat(Npt_all,Nnb_max), alpha
!
INTEGER(GI):: i, n
REAL   (DP):: q_l(Npt_all), g(Npt_all), phi(Npt_all)
REAL   (DP):: fa, fb, fc, x, y, z, a, b, c, d, e, xm, s
INTEGER(GI), PARAMETER:: Nmax = 100
REAL   (DP), PARAMETER:: E_crit = 1.D-8, E_end = 1.D-11
!
a = -1.0D-0
q_l(:) = q_i(:) + a*p(:)
CALL PHI_MESH(Npt_all, fft, psi, Qlocal, q_l, grid, phi); phi(:) = phi(:)*c_pme
CALL dUdq(Npt_all, XJZ, asym_mat, sym_mat, list, nlist, q_l, qmx, phi, g)
fa = 0.0D0; DO i=1, Npt_all; fa = fa - g(i)*p(i); END DO
b =  1.0D-0
q_l(:) = q_i(:) + b*p(:)
CALL PHI_MESH(Npt_all, fft, psi, Qlocal, q_l, grid, phi); phi(:) = phi(:)*c_pme
CALL dUdq(Npt_all, XJZ, asym_mat, sym_mat, list, nlist, q_l, qmx, phi, g)
fb = 0.0D0; DO i=1, Npt_all; fb = fb - g(i)*p(i); END DO
c = b
fc = fb
!
DO n=1, Nmax
   IF ( (fb > 0.D0 .AND. fc > 0.D0) .OR. (fb < 0.0D0 .AND. fc < 0.D0)) THEN
      c = a
      fc = fa
      d = b - a
      e = d
   END IF
   IF (DABS(fc) < DABS(fb)) THEN
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
   END IF
   xm = 0.5D0*(c - b)
   IF (DABS(xm) < E_crit .OR. DABS(fb) < E_end) THEN
      alpha = b
      !WRITE(*,100,ADVANCE = "NO") n;100 FORMAT("# Loop for Brent solver = ",I3)
      RETURN
   END IF
   IF (DABS(e) >= E_crit .AND. DABS(fa) > DABS(fb)) THEN
      s = fb/fa
      IF (a == C) THEN
         x = 2.D0*xm*s
         y = 1.D0 - s
      ELSE
         y = fa/fc
         z = fb/fc
         x = s*(2.*xm*y*(y-z) - (b-a)*(z-1.0D0))
         y = (y - 1.D0)*(z - 1.D0)*(s - 1.D0)
      END IF
      IF (x > 0.D0) y = -y
      x = DABS(x)
      IF (2.D0*x < MIN(3.*xm*y - DABS(E_crit*y), DABS(e*y))) THEN
         e = d
         d = x /y
      ELSE
         d = xm
         e = d
      END IF
   ELSE
      d = xm
      e = d
   END IF
   a = b
   fa = fb
   IF (DABS(d) > E_crit) THEN
      b = b + d
   ELSE
      b = b + SIGN(E_crit, xm)
   END IF
   q_l(:) = q_i(:) + b*p(:)
   CALL PHI_MESH(Npt_all, fft, psi, Qlocal, q_l, grid, phi); phi(:)=phi(:)*c_pme
   CALL dUdq(Npt_all, XJZ, asym_mat, sym_mat, list, nlist, q_l, qmx, phi, g)
   fb = 0.0D0; DO i=1, Npt_all; fb = fb - g(i)*p(i); END DO
END DO
WRITE(*,*) "# **** Brent solver limit breached *** "
STOP
END SUBROUTINE BRENT_SOLVER
!
!##############################################################################
SUBROUTINE dUdq(Npt_all, XJZ, asym_mat, sym_mat, list, nlist, q_i, qmx, phi, g)
USE DATAFMT
IMPLICIT NONE
INTEGER(GI):: Npt_all, list(Npt_all, Nnb_max), nlist(Npt_all)
REAL   (DP):: asym_mat(Npt_all,Nnb_max), sym_mat(Npt_all,Nnb_max), phi(Npt_all)
REAL   (DP):: XJZ(Npt_all,3), q_i(Npt_all), qmx(Npt_all,2), g(Npt_all)
!
INTEGER(GI):: i, j, id_pair
REAL   (DP):: dqmin(Npt_all), dqmax(Npt_all)
!
dqmin(:) = q_i(:) - qmx(:,1); dqmax(:) = q_i(:) - qmx(:,2)
DO i=1, Npt_all
   IF (dqmin(i) >= 0.0D0) dqmin(i) = 0.0D0
   IF (dqmax(i) <= 0.0D0) dqmax(i) = 0.0D0
END DO
g(:) = -(XJZ(:,1) + XJZ(:,2)*q_i(:) + phi(:) + 2.D0*w*( dqmin(:) + dqmax(:)) )
! XJZ(i, 2) = J_i  - 2a/\sqrt{\pi}
DO i=1, Npt_all
   DO j=1, nlist(i)
      id_pair = list(i,j)
      g(i) = g(i) - XJZ(id_pair,3)*asym_mat(i,j) - q_i(id_pair)*sym_mat(i,j)
   END DO
END DO
!
RETURN
! Returns g(:) = -dU/dq_i(:)
END SUBROUTINE dUdq
!
!##############################################################################
SUBROUTINE PHI_MESH(Npt_all, fft, psi, Qlocal, q_i, grid, phi)
USE DATAFMT
IMPLICIT NONE
!
INTEGER(GI):: Npt_all, fft(3)
REAL   (DP):: psi(fft(1),fft(2),fft(3)), Qlocal(Npt_all,4,4,4), q_i(Npt_all)
REAL   (DP):: phi(Npt_all)
INTEGER(GI):: grid(Npt_all,4,4,4,3)
!
REAL   (DP):: E_local, z(3), xfft(3)
INTEGER(GI):: i, j, k, l, ix(3)
COMPLEX(DP):: Qbas(fft(1), fft(2), fft(3)), Qinv(fft(1), fft(2), fft(3)), &
            & Qfin(fft(1), fft(2), fft(3))
!
!
xfft(:) = DBLE(fft(:))
Qbas(:,:,:) = CMPLX(0.0D0, 0.0D0)
DO i=1,Npt_all
   DO j=1,4
      DO k=1,4
         DO l=1,4
            ix(:) = grid(i,j,k,l,:)  
            Qbas(ix(1),ix(2),ix(3)) = Qbas(ix(1),ix(2),ix(3)) + &
                 & q_i(i)*CMPLX(Qlocal(i,j,k,l),0.0D0)
         END DO
      END DO
   END DO
END DO
!
Qinv(:,:,:) = Qbas(:,:,:)
CALL ZFFT3D(Qinv, fft(1), fft(2), fft(3), 0)
CALL ZFFT3D(Qinv, fft(1), fft(2), fft(3),+1)
DO i=1, fft(1)
   z(1) = (4.D0 + 2.D0*DCOS(2.D0*PI*DBLE(i-1)/xfft(1)))**2 
   DO j=1, fft(2)
      z(2) = (4.D0 + 2.D0*DCOS(2.D0*PI*DBLE(j-1)/xfft(2)))**2 
      DO k=1, fft(3)
         z(3) = (4.D0 + 2.D0*DCOS(2.D0*PI*DBLE(k-1)/xfft(3)))**2 
         Qinv(i,j,k) = Qinv(i,j,k)*PSI(i,j,k)*36.D0**3/z(1)/z(2)/z(3)
      END DO
   END DO
END DO
Qfin(:,:,:) = Qinv(:,:,:)
CALL ZFFT3D(Qfin, fft(1), fft(2), fft(3), 0)
CALL ZFFT3D(Qfin, fft(1), fft(2), fft(3),-1)
DO i=1, Npt_all
   E_local = 0.0D0
   DO j=1,4
      DO k=1,4
         DO l=1,4
            ix(:) = grid(i,j,k,l,:)
            E_local = E_local + &
                 & Qlocal(i,j,k,l)*DBLE(Qfin(ix(1),ix(2),ix(3)))
         END DO
      END DO
   END DO
   phi(i) = E_local
END DO
!
RETURN
END SUBROUTINE PHI_MESH
