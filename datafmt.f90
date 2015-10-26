MODULE DATAFMT
IMPLICIT NONE
!
! Module for data format and derived data set declaration
!
! PRECISION
! SI = single precision for integer
! DI = double precision for integer
! GI = Generic usage of integer
! SP = single precision for real
! DP = double precision for real
!http://www.lib.ncep.noaa.gov/itresources/presentations/fortran952003presentation.pdf
INTEGER, PARAMETER:: SI = SELECTED_INT_KIND(4)
INTEGER, PARAMETER:: DI = SELECTED_INT_KIND(8)
INTEGER, PARAMETER:: GI = DI 
INTEGER, PARAMETER:: SP = SELECTED_REAL_KIND(6)
INTEGER, PARAMETER:: DP = SELECTED_REAL_KIND(15)
!
! Maximum number of particles per a single cell
INTEGER(GI), PARAMETER:: Npt_max = 500
! Maximum number of neighboring particles of a single particle (in cutoff R.)
INTEGER(GI), PARAMETER:: Nnb_max = 500 
REAL   (DP), PARAMETER:: kc     = 14.399644154020082D0 
REAL   (DP), PARAMETER:: TFM    = 10.180505306645389D0
REAL   (DP), PARAMETER:: PI     = 3.1415926535897931D0
REAL   (DP), PARAMETER:: SQRTPI = 1.7724538509055159D0
REAL   (DP), PARAMETER:: kB     = 8.617343D-5
REAL   (DP), PARAMETER:: w = 20.0D0 ! unit eV/e^2 ----------- CTIP parameter
INTEGER(SI), PARAMETER:: Natom  = 5 ! Al - Ni - Fe - O - H
INTEGER(SI), PARAMETER:: Nmetal = 3 ! Al - Ni - Fe 
INTEGER(GI), PARAMETER:: Nrdf   = 200 ! Radial distribution grid
!
!
TYPE AMNT
   INTEGER(GI):: Npt_all, Npt(Natom), fft(3), Ncell_all, Ncell(3), N_o2
   INTEGER(GI):: freq_max, freq_snap, freq_rdf, freq_stat, freq_new, freq_sort
   REAL   (DP):: box(3), V, z_new, q_crit, alpha, Lcell(3)
END TYPE AMNT
!
TYPE SORT
   INTEGER(GI):: Npt, link(Npt_max), pair(13)
END type SORT
!
TYPE TCKT
   LOGICAL:: regular_run, new_o2, new_charge, snapshot, thermostat, restart
   LOGICAL:: status, cell_sort, print_corr, thermo_stoch
END TYPE TCKT
!
TYPE PTCL
   INTEGER(GI):: id
   REAL   (DP):: xx(3), xv(3), ff(3), q
END TYPE PTCL
!
TYPE POTN
   REAL(DP):: xm, qmin, qmax, chi, J, xi, Z, fafb(Natom,4)
END TYPE POTN
!
!
TYPE STAT
   REAL(DP):: box(3), Rcut, Rcut2, a, Eself, Edir, Erec, Epair, Egroup, Ectip
   REAL(DP):: T_given, T_now, mv2
   REAL(DP):: time_max, time_snap, time_rdf, time_stat
END TYPE STAT
!
TYPE METL
   REAL(DP):: r_e, f_e, rho_e, rho_s, alpha, beta, A, B, kappa, lambda, &
        & Fn(4), Fm(4), Fp(4), eta, Fe
END TYPE METL
!
TYPE OXID
   REAL(DP):: r_e(Natom), alph(Natom), beta(Natom), A(Natom), B(Natom)
   REAL(DP):: kapp(Natom), lamb(Natom), F0(4), F1(4), F2(4), F3(4), F4(4)
   REAL(DP):: re(5), rmax(5), fe, gamm, nu
END TYPE OXID
!
TYPE CORR
   INTEGER(GI):: rdf(Natom,Natom,Nrdf+1), Npt_old, Npt_new
   REAL   (DP):: dr, V
END TYPE CORR
!
CONTAINS
!
!###############################################################################
FUNCTION ADD_PTCL(x, l, n,  q)
IMPLICIT NONE
TYPE(PTCL), POINTER:: x(:), ADD_PTCL(:)
INTEGER(GI), intent(in):: l, n
TYPE(PTCL)           :: q(n)
INTEGER:: ierr
ALLOCATE(ADD_PTCL(1:l+n), STAT = ierr)
IF(ierr /=0) STOP "allocate error"
ADD_PTCL(1:l) = x(1:l)
ADD_PTCL(l+1:l+n) = q(1:n)
DEALLOCATE(x, STAT = ierr)
END FUNCTION ADD_PTCL
!
!###############################################################################
SUBROUTINE NEW_O2(Nset, sys, q)
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(STAT):: sys
TYPE(PTCL), POINTER :: q(:)
!
TYPE (PTCL):: tmp_o(Nset%N_o2*2)
REAL   (DP):: xx(Nset%N_o2*2,3), sigma, fluct, tmp, sumv(3)
INTEGER(GI):: n, i, j
!
n = Nset%N_o2*2
CALL RANDOM_POSITION_O2(n, Nset%z_new, Nset%box(1), Nset%box(2), xx)
DO i=1, n
   tmp_o(i)%xx(:) = xx(i,:)
END DO
!
WRITE(*,50)
DO i=1, n/2
   WRITE(*,100) xx(2*i-1,:), xx(2*i,:)
END DO
50 FORMAT("#")
100 FORMAT("# New oxygen dimer is spawned at (",3(F5.1,1X),") (",3(F5.1,1X),")")
sigma = 1.0D0
!tmp = DSQRT(sys%T_given/16.D0)
tmp = DSQRT(sys%T_given*3.D0/16.D0)
DO i=1,n ! Gaussian distribution of initial temperature
   DO j=1,3
      !tmp_o(i)%xv(j) = fluct(sigma)*tmp
      tmp_o(i)%xv(j) = 0.0D0
      tmp_o(i)%ff(j) = 0.0D0
   END DO
   tmp_o(i)%xv(3) = - tmp
   tmp_o(i)%q     = 0.0D0
   tmp_o(i)%id    = 4
END DO
!
q => ADD_PTCL(q, Nset%Npt_all, n, tmp_o(1:n))
Nset%Npt_all = Nset%Npt_all + n
!
RETURN 
END SUBROUTINE NEW_O2
!
!###############################################################################
SUBROUTINE POTENTIAL_SETUP(eam, param, oxygn)
IMPLICIT NONE
!
TYPE(METL):: eam(Nmetal)
TYPE(POTN):: param(Natom)
TYPE(OXID):: oxygn
!
INTEGER(GI):: i, j
REAL(DP)   :: xi, xi_a, xi_b
!
! Potential data for Embedded Atom Method concerning metal elements
!
! Al                       Ni                       Fe
eam(1)%r_e   = 2.86392D0; eam(2)%r_e   = 2.48875D0; eam(3)%r_e   = 2.48199D0
eam(1)%f_e   = 1.20378D0; eam(2)%f_e   = 2.21149D0; eam(3)%f_e   = 2.31453D0
eam(1)%rho_e =17.51747D0; eam(2)%rho_e =30.37003D0; eam(3)%rho_e =24.59573D0
eam(1)%rho_s =19.90041D0; eam(2)%rho_s =30.37137D0; eam(3)%rho_s =24.59573D0
eam(1)%alpha = 6.61317D0; eam(2)%alpha = 8.38345D0; eam(3)%alpha = 9.81827D0
eam(1)%beta  = 3.52702D0; eam(2)%beta  = 4.47117D0; eam(3)%beta  = 5.23641D0
eam(1)%A     = 0.31487D0; eam(2)%A     = 0.42905D0; eam(3)%A     = 0.39281D0
eam(1)%B     = 0.36555D0; eam(2)%B     = 0.63353D0; eam(3)%B     = 0.64624D0
eam(1)%kappa = 0.37985D0; eam(2)%kappa = 0.44360D0; eam(3)%kappa = 0.17031D0
eam(1)%lambda= 0.75969D0; eam(2)%lambda= 0.82066D0; eam(3)%lambda= 0.34061D0
eam(1)%Fn(1) =-2.80760D0; eam(2)%Fn(1) =-2.69351D0; eam(3)%Fn(1) =-2.53499D0
eam(1)%Fn(2) =-0.30144D0; eam(2)%Fn(2) =-0.07644D0; eam(3)%Fn(2) =-0.05960D0
eam(1)%Fn(3) = 1.25856D0; eam(2)%Fn(3) = 0.24144D0; eam(3)%Fn(3) = 0.19306D0
eam(1)%Fn(4) =-1.24760D0; eam(2)%Fn(4) =-2.37563D0; eam(3)%Fn(4) =-2.28232D0
eam(1)%Fm(1) =-2.83000D0; eam(2)%Fm(1) =-2.70000D0; eam(3)%Fm(1) =-2.54000D0
eam(1)%Fm(2) = 0.00000D0; eam(2)%Fm(2) = 0.00000D0; eam(3)%Fm(2) = 0.00000D0
eam(1)%Fm(3) = 0.62225D0; eam(2)%Fm(3) = 0.26539D0; eam(3)%Fm(3) = 0.20027D0
eam(1)%Fm(4) =-2.48824D0; eam(2)%Fm(4) =-0.15286D0; eam(3)%Fm(4) =-0.14877D0
eam(1)%Fp(1) =-2.83000D0; eam(2)%Fp(1) =-2.70000D0; eam(3)%Fp(1) =-2.54000D0
eam(1)%Fp(2) = 0.00000D0; eam(2)%Fp(2) = 0.00000D0; eam(3)%Fp(2) = 0.00000D0
eam(1)%Fp(3) = 0.62225D0; eam(2)%Fp(3) = 0.26539D0; eam(3)%Fp(3) = 0.20027D0
eam(1)%Fp(4) =-2.48824D0; eam(2)%Fp(4) = 4.58568D0; eam(3)%Fp(4) = 6.69465D0
eam(1)%eta   = 0.78591D0; eam(2)%eta   = 1.01318D0; eam(3)%eta   = 1.18290D0
eam(1)%Fe    =-2.82453D0; eam(2)%Fe    =-2.70839D0; eam(3)%Fe    =-2.55187D0
!
! CTIP charge data for all elements
!
! Al                       Ni                       Fe
param(1)%xm = 26.981539D0; param(2)%xm = 58.6934D0; param(3)%xm = 55.84500D0
param(1)%qmin =     0.0D0; param(2)%qmin =   0.0D0; param(3)%qmin =  0.000D0   
param(1)%qmax =     3.0D0; param(2)%qmax =   2.0D0; param(3)%qmax =  3.000D0   
param(1)%chi = -1.47914D0; param(2)%chi=-1.70804D0; param(3)%chi =-1.90587D0
param(1)%J =    9.07222D0; param(2)%J =  9.10954D0; param(3)%J =   8.99819D0
param(1)%xi =     0.968D0; param(2)%xi =   1.087D0; param(3)%xi =  1.024D000
param(1)%Z =    1.07514D0; param(2)%Z =  1.44450D0; param(3)%Z =   1.28612D0
! O                        H
param(4)%xm = 15.9994D0; param(5)%xm = 1.00794D0
param(4)%qmin =  -2.0D0; param(5)%qmin = 0.0    
param(4)%qmax =   0.0D0; param(5)%qmax = 0.0    
param(4)%chi =  2.000D0; param(5)%chi = 0.0     
param(4)%J = 14.99523D0; param(5)%J =  0.0      
param(4)%xi =   2.144D0; param(5)%xi = 0.0      
param(4)%Z =     0.00D0; param(5)%Z = 0.0       
!
! CTIP intermediate data - Coulomb integration coefficients
!
DO i=1, Natom  
   DO j=1, Natom
      IF (i == j) THEN
         xi = param(i)%xi
         param(i)%fafb(j,1) = 1.D0
         param(i)%fafb(j,2) = 11.D0*xi/8.D0
         param(i)%fafb(j,3) = 0.75D0*xi**2
         param(i)%fafb(j,4) = xi**3/6.D0
      ELSE
         xi_a = param(i)%xi
         xi_b = param(j)%xi
         param(i)%fafb(j,1) = xi_a*xi_b**4/(xi_a + xi_b)**2/(xi_a - xi_b)**2
         param(i)%fafb(j,2) = xi_b*xi_a**4/(xi_b + xi_a)**2/(xi_b - xi_a)**2 
         param(i)%fafb(j,3) = (3.D0*xi_a**2*xi_b**4 - xi_b**6) &
              & /(xi_a + xi_b)**3/(xi_a - xi_b)**3
         param(i)%fafb(j,4) = (3.D0*xi_b**2*xi_a**4 - xi_a**6) &
              & /(xi_b + xi_a)**3/(xi_b - xi_a)**3
      END IF
   END DO
END DO
!
! Oxide reaction data - with oxygen
!
! Al                       Ni                         Fe
oxygn%r_e(1)  = 2.98520D0; oxygn%r_e(2)  = 2.95732D0; oxygn%r_e(3)  = 3.07992D0
oxygn%alph(1) = 8.49741D0; oxygn%alph(2) = 7.96528D0; oxygn%alph(3) = 7.52309D0
oxygn%beta(1) = 4.52114D0; oxygn%beta(2) = 4.42411D0; oxygn%beta(3) = 4.13330D0
oxygn%A(1)    = 0.09738D0; oxygn%A(2)    = 0.13521D0; oxygn%A(3)    = 0.17108D0
oxygn%B(1)    = 0.38121D0; oxygn%B(2)    = 0.25332D0; oxygn%B(3)    = 0.39869D0
oxygn%kapp(1) = 0.18967D0; oxygn%kapp(2) = 0.47077D0; oxygn%kapp(3) = 0.22335D0
oxygn%lamb(1) = 0.95234D0; oxygn%lamb(2) = 0.65524D0; oxygn%lamb(3) = 0.34380D0
! O                        H
oxygn%r_e(4)  = 3.64857D0; oxygn%r_e(5)  = 1.0D0
oxygn%alph(4) = 5.44072D0; oxygn%alph(5) = 0.0
oxygn%beta(4) = 3.59746D0; oxygn%beta(5) = 0.0
oxygn%A(4)    = 0.34900D0; oxygn%A(5)    = 0.0
oxygn%B(4)    = 0.57438D0; oxygn%B(5)    = 0.0
oxygn%kapp(4) = 0.08007D0; oxygn%kapp(5) = 0.0
oxygn%lamb(4) = 0.39310D0; oxygn%lamb(5) = 0.0
!
oxygn%F0(1) =-1.56489D0; oxygn%F1(1) =-1.58967D0; oxygn%F2(1) =-1.54116D0
oxygn%F0(2) =-1.39123D0; oxygn%F1(2) = 1.30636D0; oxygn%F2(2) = 2.02821D0
oxygn%F0(3) = 1.77199D0; oxygn%F1(3) = 9.81033D0; oxygn%F2(3) = 6.56240D0
oxygn%F0(4) = 1.59833D0; oxygn%F1(4) = 0.00000D0; oxygn%F2(4) = 0.00000D0
!
oxygn%F3(1) =-1.51798D0; oxygn%F4(1) =-1.19082D0
oxygn%F3(2) = 2.30979D0; oxygn%F4(2) = 4.12936D0
oxygn%F3(3) = 7.69582D0; oxygn%F4(3) =10.32338D0 
oxygn%F3(4) = 0.00000D0; oxygn%F4(4) = 0.00000D0
!
oxygn%re(1) =54.62910D0; oxygn%rmax(1) =54.62910D0
oxygn%re(2) =64.26953D0; oxygn%rmax(2) =65.24078D0
oxygn%re(3) =66.21202D0; oxygn%rmax(3) =66.56797D0
oxygn%re(4) =66.92391D0; oxygn%rmax(4) =70.57748D0
oxygn%re(5) =74.23105D0; oxygn%rmax(5) =99.99999D9 ! <= just big number
!
oxygn%fe   = 1.39478D0
oxygn%gamm = 2.11725D0
oxygn%nu   = 0.37457D0
!
RETURN
END SUBROUTINE POTENTIAL_SETUP
!
END MODULE DATAFMT
