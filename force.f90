!
!#############################################################################
SUBROUTINE FORCE(Nset, eam, oxygn, sys, param, pair, q, cell)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(METL):: eam(Nmetal)
TYPE(OXID):: oxygn
TYPE(STAT):: sys
TYPE(POTN):: param(Natom)
TYPE(CORR):: pair
TYPE(PTCL):: q(Nset%Npt_all)
TYPE(SORT):: cell(Nset%Ncell_all)
!
INTEGER(GI):: i, j, m, n, l_i, l_j, id_pair, id_i, id_j, Npt_all, Ncell_all
INTEGER(GI):: Npt_local, Npt_pair, id_r
INTEGER(GI):: list(Nset%Npt_all,Nnb_max), nlist(Nset%Npt_all)
REAL   (DP):: rho(Nset%Npt_all), drhodx(Nset%Npt_all,3), fj(2), dfjdr(2)
REAL   (DP):: q_i, q_j, a, a2, box(3), phi, df, dFdrho(Nset%Npt_all)
REAL   (DP):: Epair, Egroup
REAL   (DP):: rho_s, rho_e, rho_n, rho_0, drho_list(Nset%Npt_all,Nnb_max,3)
REAL   (DP):: dx(3), rx(3), r2, rcut2, r,  DERFC, tmp, Edir, F_sum
REAL   (DP):: ifj, jfi, difjdr, djfidr, fifj, dfifjdr, xi_i, xi_j, xir_i, xir_j
REAL   (DP):: Z_i, Z_j, Ectip, exp_xir_i, exp_xir_j
!
! Initialization **************************************************************
Npt_all = Nset%Npt_all; Ncell_all = Nset%Ncell_all
DO i=1, Npt_all; q(i)%ff(:) = 0.0D0; END DO
rho(:) = 0.0D0
drhodx(:,:) = 0.0D0
box(:) = sys%box(:)
rcut2 = sys%Rcut2
a = sys%a
a2 = a**2
Edir = 0.0D0; Epair = 0.0D0; Egroup = 0.0D0; Ectip = 0.0D0
nlist(:) = 0
list(:,:) = 0
drho_list(:,:,:) = 0.0D0
pair%Npt_new = Npt_all
!
DO m=1, Ncell_all
   !
   ! Interactions inside of a single cell
   Npt_local = cell(m)%Npt
   DO i=1, Npt_local-1
      l_i = cell(m)%link(i)
      id_i = q(l_i)%id;      q_i = q(l_i)%q
      xi_i = param(id_i)%xi; Z_i = param(id_i)%Z
      DO j=i+1, Npt_local
         l_j = cell(m)%link(j)
         dx(:) = q(l_i)%xx(:) - q(l_j)%xx(:)
         rx(:) = dx(:) - box(:)*DNINT(dx(:)/box(:))
         r2 = rx(1)**2 + rx(2)**2 + rx(3)**2
         IF (r2 < rcut2) THEN
            r = DSQRT(r2); id_r = INT(r/pair%dr)
            id_j = q(l_j)%id;      q_j = q(l_j)%q
            xi_j = param(id_j)%xi; Z_j = param(id_j)%Z
            IF ((id_i < 4 .AND. id_j < 4) .AND. id_i == id_j) THEN
               CALL METAL_HOMO(eam(id_i), r, fj, dfjdr, phi, df)
            ELSE IF (id_i==1 .AND. id_j==2) THEN
               CALL AL_NI(eam(id_i), eam(id_j), r, fj, dfjdr, phi, df)
            ELSE IF (id_i==2 .AND. id_j==1) THEN
               CALL AL_NI(eam(id_j), eam(id_i), r, fj, dfjdr, phi, df)
            ELSE IF ((id_i < 4 .AND. id_j < 4) .AND. id_i /= id_j) THEN
               CALL METAL_HETERO(eam(id_i), eam(id_j), r, fj, dfjdr, phi, df)
            ELSE IF (id_i == 4 .AND. id_j < 4) THEN
               CALL OXYGEN_METAL(oxygn, eam(id_j), id_j, r, fj, dfjdr, phi, df)
            ELSE IF (id_i < 4 .AND. id_j == 4) THEN
               CALL METAL_OXYGEN(eam(id_i), id_i, oxygn, r, fj, dfjdr, phi, df)
            ELSE IF (id_i == 4 .AND. id_j == 4) THEN
               CALL OXYGEN_ONLY(oxygn, r, fj, dfjdr, phi, df)
            ELSE
               STOP "=== Unknown particle pair ==="
            END IF
            !
            pair%rdf(id_i,id_j,id_r) = pair%rdf(id_i,id_j,id_r) + 2
            nlist(l_i) = nlist(l_i) + 1; nlist(l_j) = nlist(l_j) + 1
            IF (nlist(i) > Nnb_max .OR. nlist(j) > Nnb_max) &
                 & STOP "=== Memeory error for EAM density ==="
            list(l_i,nlist(l_i)) = l_j; list(l_j,nlist(l_j)) = l_i
            drho_list(l_i,nlist(l_i),:) =  rx(:)*dfjdr(1)
            drho_list(l_j,nlist(l_j),:) = -rx(:)*dfjdr(2)
            !
            rho(l_i) = rho(l_i) + fj(2)
            rho(l_j) = rho(l_j) + fj(1)
            drhodx(l_i,:) = drhodx(l_i,:) + rx(:)*dfjdr(2)
            drhodx(l_j,:) = drhodx(l_j,:) - rx(:)*dfjdr(1)
            !
            q(l_i)%ff(:) = q(l_i)%ff(:) + rx(:)*df
            q(l_j)%ff(:) = q(l_j)%ff(:) - rx(:)*df
            Epair = Epair + phi
            !
            ! Direct summation for Ewald (PME) ********************************
            tmp = DERFC(a*r)
            Edir = Edir + q_i * q_j * tmp/r
            df = kc * q_i * q_j * (2.D0*DEXP(-a2*r2)*a*r/SQRTPI + tmp)/r2/r
            q(l_i)%ff(:) = q(l_i)%ff(:) + df*rx(:)
            q(l_j)%ff(:) = q(l_j)%ff(:) - df*rx(:)
            !
            ! Other electrostatic components from CTIP ************************
            IF (id_i == id_j) THEN
               xir_i   = 2.0*xi_i*r;  exp_xir_i = DEXP(-xir_i)
               jfi     = -exp_xir_i*(xi_i + 1.D0/r)
               ifj     = jfi
               djfidr  = exp_xir_i*(2.D0*xi_i**2 + (1.D0 + xir_i)/r2)
               difjdr  = djfidr
               fifj    = - exp_xir_i*(1.D0 + &
                    &     param(id_i)%fafb(id_j,2)*r + &
                    &     param(id_i)%fafb(id_j,3)*r2 + &
                    &     param(id_i)%fafb(id_j,4)*r2*r)
               dfifjdr = exp_xir_i * ( &
                    & (xir_i + 1.D0)*(1.D0/r2 + &
                    &                  param(id_i)%fafb(id_j,2)/r + &
                    &                  param(id_i)%fafb(id_j,3) + &
                    &                  param(id_i)%fafb(id_j,4)*r) &
                    & - ( param(id_i)%fafb(id_j,2)/r + &
                    &     param(id_i)%fafb(id_j,3)*2.D0 + &
                    &     param(id_i)%fafb(id_j,4)*3.D0*r))
            ELSE
               xir_i   = 2.D0*xi_i*r; exp_xir_i = DEXP(-xir_i)
               xir_j   = 2.D0*xi_j*r; exp_xir_j = DEXP(-xir_j)
               jfi     = -exp_xir_i*(xi_i + 1.D0/r)
               ifj     = -exp_xir_j*(xi_j + 1.D0/r)
               djfidr  = exp_xir_i*(2.D0*xi_i**2 + (1.D0 + xir_i)/r2)
               difjdr  = exp_xir_j*(2.D0*xi_j**2 + (1.D0 + xir_j)/r2)
               fifj    = - param(id_i)%fafb(id_j,1)*exp_xir_i - &
                    &      param(id_i)%fafb(id_j,2)*exp_xir_j - &
                    &      param(id_i)%fafb(id_j,3)*exp_xir_i/r - &
                    &      param(id_i)%fafb(id_j,4)*exp_xir_j/r
               dfifjdr = param(id_i)%fafb(id_j,1)*exp_xir_i*2.D0*xi_i + &
                    &    param(id_i)%fafb(id_j,2)*exp_xir_j*2.D0*xi_j + &
                    &    param(id_i)%fafb(id_j,3)*exp_xir_i*(1.D0 + xir_i)/r2+ &
                    &    param(id_i)%fafb(id_j,4)*exp_xir_j*(1.D0 + xir_j)/r2
            END IF
            df = kc*( (q_i*Z_j + q_j*Z_i - q_i*q_j) * dfifjdr - &
                 & q_i*Z_j*djfidr - q_j*Z_i*difjdr )/r
            q(l_i)%ff(:) = q(l_i)%ff(:) + df*rx(:)
            q(l_j)%ff(:) = q(l_j)%ff(:) - df*rx(:)
            Ectip = Ectip + q_i*Z_j*(jfi - fifj) + q_j*Z_i*(ifj - fifj) + &
                 &  q_i*q_j*fifj
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
         id_j = q(l_j)%id;      q_j = q(l_j)%q
         xi_j = param(id_j)%xi; Z_j = param(id_j)%Z
         DO i=1, Npt_local
            l_i = cell(m)%link(i)
            dx(:) = q(l_i)%xx(:) - q(l_j)%xx(:)
            rx(:) = dx(:) - box(:)*DNINT(dx(:)/box(:))
            r2 = rx(1)**2 + rx(2)**2 + rx(3)**2
            IF (r2 < rcut2) THEN
               r = DSQRT(r2); id_r = INT(r/pair%dr)
               id_i = q(l_i)%id;      q_i = q(l_i)%q
               xi_i = param(id_i)%xi; Z_i = param(id_i)%Z
               IF ((id_i < 4 .AND. id_j < 4) .AND. id_i == id_j) THEN
                  CALL METAL_HOMO(eam(id_i), r, fj, dfjdr, phi, df)
               ELSE IF (id_i==1 .AND. id_j==2) THEN
                  CALL AL_NI(eam(id_i), eam(id_j), r, fj, dfjdr, phi, df)
               ELSE IF (id_i==2 .AND. id_j==1) THEN
                  CALL AL_NI(eam(id_j), eam(id_i), r, fj, dfjdr, phi, df)
               ELSE IF ((id_i < 4 .AND. id_j < 4) .AND. id_i /= id_j) THEN
                  CALL METAL_HETERO(eam(id_i), eam(id_j), r, fj, dfjdr, phi, df)
               ELSE IF (id_i == 4 .AND. id_j < 4) THEN
                  CALL OXYGEN_METAL(oxygn, eam(id_j), id_j, r, fj, dfjdr,phi,df)
               ELSE IF (id_i < 4 .AND. id_j == 4) THEN
                  CALL METAL_OXYGEN(eam(id_i), id_i, oxygn, r, fj, dfjdr,phi,df)
               ELSE IF (id_i == 4 .AND. id_j == 4) THEN
                  CALL OXYGEN_ONLY(oxygn, r, fj, dfjdr, phi, df)
               ELSE
                  STOP "=== Unknown particle pair ==="
               END IF
               !
               pair%rdf(id_i,id_j,id_r) = pair%rdf(id_i,id_j,id_r) + 2
               nlist(l_i) = nlist(l_i) + 1; nlist(l_j) = nlist(l_j) + 1
               IF (nlist(i) > Nnb_max .OR. nlist(j) > Nnb_max) &
                    & STOP "=== Memeory error for EAM density ==="
               list(l_i,nlist(l_i)) = l_j; list(l_j,nlist(l_j)) = l_i
               drho_list(l_i,nlist(l_i),:) =  rx(:)*dfjdr(1)
               drho_list(l_j,nlist(l_j),:) = -rx(:)*dfjdr(2)
               !
               rho(l_i) = rho(l_i) + fj(2)
               rho(l_j) = rho(l_j) + fj(1)
               drhodx(l_i,:) = drhodx(l_i,:) + rx(:)*dfjdr(2)
               drhodx(l_j,:) = drhodx(l_j,:) - rx(:)*dfjdr(1)
               !
               q(l_i)%ff(:) = q(l_i)%ff(:) + rx(:)*df
               q(l_j)%ff(:) = q(l_j)%ff(:) - rx(:)*df
               Epair = Epair + phi
               !
               ! Direct summation for Ewald (PME) *****************************
               tmp = DERFC(a*r)
               Edir = Edir + q_i * q_j * tmp/r
               df = kc * q_i * q_j * (2.D0*DEXP(-a2*r2)*a*r/SQRTPI + tmp)/r2/r
               q(l_i)%ff(:) = q(l_i)%ff(:) + df*rx(:)
               q(l_j)%ff(:) = q(l_j)%ff(:) - df*rx(:)
               !
               ! Other electrostatic components from CTIP *********************
               IF (id_i == id_j) THEN
                  xir_i   = 2.0*xi_i*r;  exp_xir_i = DEXP(-xir_i)
                  jfi     = -exp_xir_i*(xi_i + 1.D0/r)
                  ifj     = jfi
                  djfidr  = exp_xir_i*(2.D0*xi_i**2 + (1.D0 + xir_i)/r2)
                  difjdr  = djfidr
                  fifj    = - exp_xir_i*(1.D0 + &
                       &     param(id_i)%fafb(id_j,2)*r + &
                       &     param(id_i)%fafb(id_j,3)*r2 + &
                       &     param(id_i)%fafb(id_j,4)*r2*r)
                  dfifjdr = exp_xir_i * ( &
                       & (xir_i + 1.D0)*(1.D0/r2 + &
                       &                  param(id_i)%fafb(id_j,2)/r + &
                       &                  param(id_i)%fafb(id_j,3) + &
                       &                  param(id_i)%fafb(id_j,4)*r) &
                       & - ( param(id_i)%fafb(id_j,2)/r + &
                       &     param(id_i)%fafb(id_j,3)*2.D0 + &
                       &     param(id_i)%fafb(id_j,4)*3.D0*r))
               ELSE
                  xir_i   = 2.D0*xi_i*r; exp_xir_i = DEXP(-xir_i)
                  xir_j   = 2.D0*xi_j*r; exp_xir_j = DEXP(-xir_j)
                  jfi     = -exp_xir_i*(xi_i + 1.D0/r)
                  ifj     = -exp_xir_j*(xi_j + 1.D0/r)
                  djfidr  = exp_xir_i*(2.D0*xi_i**2 + (1.D0 + xir_i)/r2)
                  difjdr  = exp_xir_j*(2.D0*xi_j**2 + (1.D0 + xir_j)/r2)
                  fifj    = - param(id_i)%fafb(id_j,1)*exp_xir_i - &
                       &      param(id_i)%fafb(id_j,2)*exp_xir_j - &
                       &      param(id_i)%fafb(id_j,3)*exp_xir_i/r - &
                       &      param(id_i)%fafb(id_j,4)*exp_xir_j/r
                  dfifjdr = param(id_i)%fafb(id_j,1)*exp_xir_i*2.D0*xi_i + &
                       &    param(id_i)%fafb(id_j,2)*exp_xir_j*2.D0*xi_j + &
                       &    param(id_i)%fafb(id_j,3)*exp_xir_i*(1.D0+xir_i)/r2+&
                       &    param(id_i)%fafb(id_j,4)*exp_xir_j*(1.D0+xir_j)/r2
               END IF
               df = kc*( (q_i*Z_j + q_j*Z_i - q_i*q_j) * dfifjdr - &
                    & q_i*Z_j*djfidr - q_j*Z_i*difjdr )/r
               q(l_i)%ff(:) = q(l_i)%ff(:) + df*rx(:)
               q(l_j)%ff(:) = q(l_j)%ff(:) - df*rx(:)
               Ectip = Ectip + q_i*Z_j*(jfi - fifj) + q_j*Z_i*(ifj - fifj) + &
                    &  q_i*q_j*fifj
            END IF
         END DO
      END DO
   END DO
END DO
sys%Epair = Epair
sys%Edir = Edir * kc
sys%Ectip = Ectip * kc
dFdrho(:) = 0.0D0
!
!
DO i=1, Npt_all
   id_i = q(i)%id
   !
   ! Metal only
   IF (id_i < 4) THEN
      rho_s = eam(id_i)%rho_s
      rho_e = eam(id_i)%rho_e      
      rho_n = 0.85D0*rho_e
      rho_0 = 1.15D0*rho_e
      F_sum = 0.0D0
      IF (rho(i) < rho_n) THEN
         tmp = rho(i)/rho_n - 1.D0
         DO j=1,4
            F_sum = F_sum + eam(id_i)%Fn(j)*tmp**(j-1)
         END DO
         DO j=2,4
            dFdrho(i) = dFdrho(i) + eam(id_i)%Fn(j)*DBLE(j-1)*tmp**(j-2)/rho_n
         END DO
      ELSE IF (rho(i) < rho_e) THEN
         tmp = rho(i)/rho_e - 1.D0
         DO j=1,4
            F_sum = F_sum + eam(id_i)%Fm(j)*tmp**(j-1)
         END DO
         DO j=2,4
            dFdrho(i) = dFdrho(i) + eam(id_i)%Fm(j)*DBLE(j-1)*tmp**(j-2)/rho_e
         END DO
      ELSE IF (rho(i) < rho_0) THEN
         tmp = rho(i)/rho_e - 1.D0
         DO j=1,4
            F_sum = F_sum + eam(id_i)%Fp(j)*tmp**(j-1)
         END DO
         DO j=2,4
            dFdrho(i) = dFdrho(i) + eam(id_i)%Fp(j)*DBLE(j-1)*tmp**(j-2)/rho_e
         END DO
      ELSE IF (rho(i) > 1.D-6) THEN
         tmp = (rho(i)/rho_s)**eam(id_i)%eta
         F_sum = eam(id_i)%Fe*(1.D0 - DLOG(tmp))*tmp
         dFdrho(i) = - eam(id_i)%Fe*eam(id_i)%eta*DLOG(tmp)*tmp/rho(i)
      END IF
      !
      ! Oxygen only
   ELSE IF (id_i == 4) THEN
      F_sum = 0.0D0
      IF (rho(i) < oxygn%rmax(1)) THEN
         rho_e = oxygn%re(1)
         tmp = rho(i)/rho_e - 1.D0
         DO j=1,4
            F_sum = F_sum + oxygn%F0(j)*tmp**(j-1)
         END DO
         DO j=2,4
            dFdrho(i) = dFdrho(i) + oxygn%F0(j)*DBLE(j-1)*tmp**(j-2)/rho_e
         END DO
      ELSE IF (rho(i) < oxygn%rmax(2)) THEN
         rho_e = oxygn%re(2)
         tmp = rho(i)/rho_e - 1.D0
         DO j=1,4
            F_sum = F_sum + oxygn%F1(j)*tmp**(j-1)
         END DO
         DO j=2,4
            dFdrho(i) = dFdrho(i) + oxygn%F1(j)*DBLE(j-1)*tmp**(j-2)/rho_e
         END DO
      ELSE IF (rho(i) < oxygn%rmax(3)) THEN
         rho_e = oxygn%re(3)
         tmp = rho(i)/rho_e - 1.D0
         DO j=1,4
            F_sum = F_sum + oxygn%F2(j)*tmp**(j-1)
         END DO
         DO j=2,4
            dFdrho(i) = dFdrho(i) + oxygn%F2(j)*DBLE(j-1)*tmp**(j-2)/rho_e
         END DO
      ELSE IF (rho(i) < oxygn%rmax(4)) THEN
         rho_e = oxygn%re(4)
         tmp = rho(i)/rho_e - 1.D0
         DO j=1,4
            F_sum = F_sum + oxygn%F3(j)*tmp**(j-1)
         END DO
         DO j=2,4
            dFdrho(i) = dFdrho(i) + oxygn%F3(j)*DBLE(j-1)*tmp**(j-2)/rho_e
         END DO
      ELSE  
         rho_e = oxygn%re(5)
         tmp = rho(i)/rho_e - 1.D0
         DO j=1,4
            F_sum = F_sum + oxygn%F4(j)*tmp**(j-1)
         END DO
         DO j=2,4
            dFdrho(i) = dFdrho(i) + oxygn%F4(j)*DBLE(j-1)*tmp**(j-2)/rho_e
         END DO
      END IF
   ELSE 
      PRINT *, "# Unknown element", id_i
      STOP
   END IF
   q(i)%ff(:) = q(i)%ff(:) - drhodx(i,:)*dFdrho(i)
   Egroup = Egroup + F_sum
END DO
sys%Egroup = Egroup
DO i=1, Npt_all
   DO j = 1, nlist(i)
      id_j = list(i,j)
      q(i)%ff(:) = q(i)%ff(:) - dFdrho(id_j)*drho_list(i,j,:)
   END DO
END DO
!
RETURN
END SUBROUTINE FORCE
!
!##############################################################################
SUBROUTINE PME(Nset, sys, q, psi)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(STAT):: sys
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: psi(Nset%fft(1),Nset%fft(2),Nset%fft(3))
!
INTEGER(GI):: i, j, k, l, id, fft(3), nn(3), ix(3), Npt_all
REAL   (DP):: V, q_i, ff(Nset%Npt_all,3), u(3), z(3), box(3), xfft(3), kx(3)
REAL   (DP):: AX(Nset%Npt_all,4,3), AXX(Nset%Npt_all,4,3), dQdx(3), const, Epme
REAL   (DP):: Fsum(3)
COMPLEX(DP):: Qbas(Nset%fft(1), Nset%fft(2), Nset%fft(3)), &
     & Qinv(Nset%fft(1), Nset%fft(2), Nset%fft(3)), &
     & Qfin(Nset%fft(1), Nset%fft(2), Nset%fft(3))
!
Npt_all = Nset%Npt_all
V = Nset%box(1)*Nset%box(2)*Nset%box(3)
box(:) = Nset%box(:); xfft(:) = DBLE(Nset%fft(:)); fft(:) = Nset%fft(:)
Qbas(:,:,:) = CMPLX(0.0D0, 0.0D0)
ff(:,:) = 0.0D0
DO i=1, Npt_all
   id = q(i)%id
   q_i = q(i)%q
   u(:) = xfft(:)*(q(i)%xx(:) + box(:)*0.5D0)/box(:)
   z(:) = DINT(u(:)) - u(:) + 1.D0
   AX(i,1,:)  =       z(:)**3                                      /6.D0
   AXX(i,1,:) = -3.D0*Z(:)**2                                      /6.D0
   z(:) = z(:) + 1.D0
   AX(i,2,:) = (-3.D0*z(:)**3 + 12.D0*z(:)**2 - 12.D0*z(:) +  4.D0)/6.D0
   AXX(i,2,:)= ( 9.D0*z(:)**2 - 24.D0*z(:)    + 12.D0             )/6.D0
   z(:) = z(:) + 1.D0
   AX(i,3,:) = ( 3.D0*z(:)**3 - 24.D0*z(:)**2 + 60.D0*z(:) - 44.D0)/6.D0
   AXX(i,3,:)= (-9.D0*z(:)**2 + 48.D0*z(:)    - 60.D0             )/6.D0
   z(:) = z(:) + 1.D0
   AX(i,4,:) = (     -z(:)**3 + 12.D0*z(:)**2 - 48.D0*z(:) + 64.D0)/6.D0
   AXX(i,4,:)= ( 3.D0*z(:)**2 - 24.D0*z(:)    + 48.D0             )/6.D0
   DO j=1,4; nn(1) = j
      DO k=1,4; nn(2) = k
         DO l=1,4; nn(3) = l
            ix(:) = MOD(INT(u(:)) + nn(:) - 2 + fft(:), fft(:)) + 1
            Qbas(ix(1),ix(2),ix(3)) = Qbas(ix(1),ix(2),ix(3)) + &
                 & q_i*CMPLX(AX(i,nn(1),1)*AX(i,nn(2),2)*AX(i,nn(3),3),0.0D0)
         END DO
      END DO
   END DO
END DO
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
Epme = 0.0D0
DO i=1, fft(1)
   DO j=1, fft(2)
      DO k=1, fft(3)
         Epme = Epme + Qbas(i,j,k)*Qfin(i,j,k)
      END DO
   END DO
END DO
Epme = Epme * xfft(1) * xfft(2) * xfft(3)
const = kc * xfft(1) * xfft(2)* xfft(3) / V / PI
kx(:) = xfft(:)/box(:)
Fsum(:) = 0.0D0
DO i=1, Npt_all
   id = q(i)%id
   q_i = q(i)%q
   u(:) = xfft(:)*(q(i)%xx(:) + box(:)*0.5D0)/box(:)
   DO j=1,4; nn(1) = j
      DO k=1,4; nn(2) = k
         DO l=1,4; nn(3) = l
            dQdx(1) = q_i*AXX(i,nn(1),1)* AX(i,nn(2),2)* AX(i,nn(3),3)*kx(1)
            dQdx(2) = q_i* AX(i,nn(1),1)*AXX(i,nn(2),2)* AX(i,nn(3),3)*kx(2)
            dQdx(3) = q_i* AX(i,nn(1),1)* AX(i,nn(2),2)*AXX(i,nn(3),3)*kx(3)

            ix(:) = MOD(INT(u(:)) + nn(:) - 2 + fft(:), fft(:)) + 1
            ff(i,:) = ff(i,:) - const*dQdx(:)*DBLE(Qfin(ix(1),ix(2),ix(3)))
         END DO
      END DO
   END DO
   Fsum(:) = Fsum(:) + ff(i,:)
END DO
Fsum(:) = Fsum(:)/DBLE(Npt_all)
DO i=1, Npt_all
   q(i)%ff(:) = q(i)%ff(:) + ff(i,:) - Fsum(:)
END DO
sys%Erec = kc * Epme / V / PI / 2.D0
!
RETURN
END SUBROUTINE PME
!
!##############################################################################
SUBROUTINE PME_SELF(Nset, sys, param, q, psi)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(STAT):: sys
TYPE(POTN):: param(Natom)
TYPE(PTCL):: q(Nset%Npt_all)
REAL(DP)  :: psi(Nset%fft(1), Nset%fft(2), Nset%fft(3))

!
REAL   (DP):: m2, q_i, Eself, a, a2, X, Jx, dq, Esum
INTEGER(GI):: i, j, k, ix, iy, iz, m, Npt_all, id
!
Npt_all = Nset%Npt_all
a = 3.8D0/sys%Rcut
sys%a = a
a2 = sys%a**2
Eself = 0.0D0
Esum = 0.0D0
DO i=1, Npt_all
   q_i = q(i)%q
   Eself = Eself - q_i**2
   id = q(i)%id
   X  = param(id)%chi
   Jx = param(id)%J
   Esum = Esum + X*q_i + 0.5D0*Jx*q_i**2
   dq = q_i - param(id)%qmin
   IF (dq < 0.0D0) Esum = Esum + 2.D0*w*dq**2
   dq = q_i - param(id)%qmax
   IF (dq > 0.0D0) Esum = Esum + 2.D0*w*dq**2   
END DO
sys%Eself = Eself * a / SQRTPI * kc + Esum
!
DO i=1, Nset%fft(1)
   ix = MOD(i - 1 + Nset%fft(1)/2, Nset%fft(1)) - Nset%fft(1)/2
   DO j=1, Nset%fft(2)
      iy = MOD(j - 1 + Nset%fft(2)/2, Nset%fft(2)) - Nset%fft(2)/2
      DO k=1, Nset%fft(3)
         iz = MOD(k - 1 + Nset%fft(3)/2, Nset%fft(3)) - Nset%fft(3)/2
         m = ix**2 + iy**2 + iz**2
         IF (m /= 0) THEN
            m2 = DBLE(ix**2) / sys%box(1)**2 + DBLE(iy**2) / sys%box(2)**2 + &
                 & DBLE(iz**2) / sys%box(3)**2
            psi(i,j,k) = DEXP(-m2*PI**2/a2)/m2
         END IF
      END DO
   END DO
END DO
!
RETURN
END SUBROUTINE PME_SELF
         
