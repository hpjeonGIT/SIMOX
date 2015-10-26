!
!#############################################################################
SUBROUTINE METAL_HOMO(eam, r, fj, dfjdr, phi, df)
USE DATAFMT
IMPLICIT NONE
TYPE(METL):: eam
REAL  (DP):: r, fj(2), dfjdr(2), phi, df
!
REAL   (DP):: re, fe, alpha, beta, lambda, kappa, A, B
REAL   (DP):: r_re, r_re_1, ftn, dftndr, dphidr
REAL   (DP):: lambda20, lambda19, kappa20, kappa19, exp_alpha, exp_beta
!
re     = eam%r_e;   fe    = eam%f_e;   alpha  = eam%alpha; beta  = eam%beta; 
lambda = eam%lambda;kappa = eam%kappa; A      = eam%A;     B     = eam%B
!
r_re = r/re; r_re_1 = r_re - 1.D0
lambda20 = (r_re - lambda)**20; lambda19 = lambda20/(r_re - lambda)
kappa20  = (r_re - kappa)**20;  kappa19  = kappa20/(r_re - kappa)
exp_alpha = DEXP(-alpha*r_re_1)
exp_beta  = DEXP(-beta*r_re_1)
ftn  = exp_beta / (1.D0 + lambda20)
dftndr = exp_beta*(-beta*(1.D0 + lambda20) - 20.D0*lambda19) &
     &   /re/(1.D0 + lambda20)**2
!
fj(1) = fe*ftn
fj(2) = fj(1)
dfjdr(1) = fe*dftndr/r
dfjdr(2) = dfjdr(1)
phi = A*exp_alpha/(1.D0 + kappa20) - B*ftn
dphidr = A*exp_alpha*(-alpha*(1.D0 + kappa20) - 20.D0*kappa19) &
     &   /re/(1.D0 + kappa20)**2 - B*dftndr
df = -dphidr/r
!
RETURN
END SUBROUTINE METAL_HOMO
!
!#############################################################################
SUBROUTINE METAL_HETERO(eam_i, eam_j, r, fj, dfjdr, phi, df)
USE DATAFMT
IMPLICIT NONE
TYPE(METL):: eam_i, eam_j
REAL  (DP):: r, fj(2), dfjdr(2), phi, df
!
REAL   (DP):: re, fe, alpha, beta, lambda, kappa, A, B
REAL   (DP):: r_re, r_re_1, ftn, dftndr, phi_i, phi_j, dphi_i, dphi_j
REAL   (DP):: lambda20, lambda19, kappa20, kappa19, exp_alpha, exp_beta
!
re     = eam_i%r_e;    fe    = eam_i%f_e; 
alpha  = eam_i%alpha;  beta  = eam_i%beta; 
lambda = eam_i%lambda; kappa = eam_i%kappa; 
A      = eam_i%A;      B     = eam_i%B
!
r_re = r/re; r_re_1 = r_re - 1.D0
lambda20  = (r_re - lambda)**20; lambda19 = lambda20/(r_re - lambda)
kappa20   = (r_re - kappa)**20;  kappa19  = kappa20/(r_re - kappa)
exp_alpha = DEXP(-alpha*r_re_1)
exp_beta  = DEXP(-beta*r_re_1)
ftn       = exp_beta / (1.D0 + lambda20)
dftndr    = exp_beta*(-beta*(1.D0 + lambda20) - 20.D0*lambda19) &
     &      /re/(1.D0 + lambda20)**2
!
! *************************** part_i
fj(1) = fe*ftn
dfjdr(1) = fe*dftndr/r
phi_i = A*exp_alpha/(1.D0 + kappa20) - B*ftn
dphi_i = A*exp_alpha*(-alpha*(1.D0 + kappa20) - 20.D0*kappa19) &
     &   /re/(1.D0 + kappa20)**2 - B*dftndr
!
re     = eam_j%r_e;   fe    = eam_j%f_e;   
alpha  = eam_j%alpha; beta  = eam_j%beta; 
lambda = eam_j%lambda;kappa = eam_j%kappa; 
A      = eam_j%A;     B     = eam_j%B
r_re = r/re; r_re_1 = r_re - 1.D0
lambda20  = (r_re - lambda)**20; lambda19 = lambda20/(r_re - lambda)
kappa20   = (r_re - kappa)**20;  kappa19  = kappa20/(r_re - kappa)
exp_alpha = DEXP(-alpha*r_re_1)
exp_beta  = DEXP(-beta*r_re_1)
ftn       = exp_beta / (1.D0 + lambda20)
dftndr    = exp_beta*(-beta*(1.D0 + lambda20) - 20.D0*lambda19) &
     &      /re/(1.D0 + lambda20)**2
!
! **************************** part j
fj(2)    = fe*ftn
dfjdr(2) = fe*dftndr/r
phi_j = A*exp_alpha/(1.D0 + kappa20) - B*ftn
dphi_j = A*exp_alpha*(-alpha*(1.D0 + kappa20) - 20.D0*kappa19) &
     &   /re/(1.D0 + kappa20)**2 - B*dftndr
!
! **************************** hybrid
phi = 0.5D0*(phi_i*fj(2)/fj(1) + phi_j*fj(1)/fj(2))
df = - 0.5D0*( &
     & (dphi_i*fj(2)*fj(1) + phi_i*dfjdr(2)*r*fj(1) - &
     &  phi_i*fj(2)*dfjdr(1)*r)/fj(1)**2 + &
     & (dphi_j*fj(1)*fj(2) + phi_j*dfjdr(1)*r*fj(2) - &
     &  phi_j*fj(1)*dfjdr(2)*r)/fj(2)**2 )/r
!
RETURN
END SUBROUTINE METAL_HETERO
!
!#############################################################################
SUBROUTINE AL_NI(eam_i, eam_j, r, fj, dfjdr, phi, df)
USE DATAFMT
IMPLICIT NONE
TYPE(METL):: eam_i, eam_j
REAL  (DP):: r, fj(2), dfjdr(2), phi, df
!
REAL   (DP):: re, fe, alpha, beta, lambda, kappa, A, B
REAL   (DP):: r_re, r_re_1, ftn, dftndr, dphidr
REAL   (DP):: lambda20, lambda19, kappa20, kappa19, exp_alpha, exp_beta
!
re     = eam_i%r_e;    
fe     = eam_i%f_e; 
beta   = eam_i%beta; 
lambda = eam_i%lambda;
!
r_re = r/re; r_re_1 = r_re - 1.D0
lambda20  = (r_re - lambda)**20; lambda19 = lambda20/(r_re - lambda)
kappa20   = (r_re - kappa)**20;  kappa19  = kappa20/(r_re - kappa)
exp_alpha = DEXP(-alpha*r_re_1)
exp_beta  = DEXP(-beta*r_re_1)
ftn       = exp_beta / (1.D0 + lambda20)
dftndr    = exp_beta*(-beta*(1.D0 + lambda20) - 20.D0*lambda19) &
     &      /re/(1.D0 + lambda20)**2
!
! *************************** part_i
fj(1) = fe*ftn
dfjdr(1) = fe*dftndr/r
!
!
re     = eam_j%r_e;   
fe     = eam_j%f_e;   
beta   = eam_j%beta; 
lambda = eam_j%lambda;
!
r_re = r/re; r_re_1 = r_re - 1.D0
lambda20  = (r_re - lambda)**20; lambda19 = lambda20/(r_re - lambda)
kappa20   = (r_re - kappa)**20;  kappa19  = kappa20/(r_re - kappa)
exp_alpha = DEXP(-alpha*r_re_1)
exp_beta  = DEXP(-beta*r_re_1)
ftn       = exp_beta / (1.D0 + lambda20)
dftndr    = exp_beta*(-beta*(1.D0 + lambda20) - 20.D0*lambda19) &
     &      /re/(1.D0 + lambda20)**2
!
! **************************** part j
fj(2)    = fe*ftn
dfjdr(2) = fe*dftndr/r
!
re     = 2.71579D0
alpha  = 8.00443D0
beta   = 4.75970D0
lambda = 0.81777D0
kappa  = 0.63279D0
A      = 0.44254D0
B      = 0.68349D0
r_re = r/re; r_re_1 = r_re - 1.D0
lambda20  = (r_re - lambda)**20; lambda19 = lambda20/(r_re - lambda)
kappa20   = (r_re - kappa)**20;  kappa19  = kappa20/(r_re - kappa)
exp_alpha = DEXP(-alpha*r_re_1)
exp_beta  = DEXP(-beta*r_re_1)
phi = A*exp_alpha/(1.D0 + kappa20) - B*exp_beta / (1.D0 + lambda20)
dphidr = A*exp_alpha*(-alpha*(1.D0 + kappa20) - 20.D0*kappa19) &
     &   /re/(1.D0 + kappa20)**2 - &
     &   B*exp_beta*(-beta*(1.D0 + lambda20) - 20.D0*lambda19) &
     &   /re/(1.D0 + lambda20)**2
df = -dphidr/r
!
RETURN
END SUBROUTINE AL_NI
!
!#############################################################################
SUBROUTINE METAL_OXYGEN(eam, id_m, oxygn, r, fj, dfjdr, phi, df)
USE DATAFMT
IMPLICIT NONE
TYPE(METL) :: eam
TYPE(OXID) :: oxygn
REAL  (DP) :: r, fj(2), dfjdr(2), phi, df
INTEGER(GI):: id_m
!
REAL   (DP):: re, fe, alpha, beta, lambda, kappa, A, B
REAL   (DP):: r_re, r_re_1, ftn, dftndr, dphidr, gamm, nu, nu20, nu19
REAL   (DP):: lambda20, lambda19, kappa20, kappa19, exp_alpha, exp_beta, exp_g
!
re     = eam%r_e;   fe    = eam%f_e;   beta  = eam%beta; lambda = eam%lambda
!
r_re = r/re; r_re_1 = r_re - 1.D0
lambda20 = (r_re - lambda)**20; lambda19 = lambda20/(r_re - lambda)
exp_beta  = DEXP(-beta*r_re_1)
ftn  = exp_beta / (1.D0 + lambda20)
dftndr = exp_beta*(-beta*(1.D0 + lambda20) - 20.D0*lambda19) &
     &   /re/(1.D0 + lambda20)**2
!
! ********************** METAL part (electron density only)
fj(1) = fe*ftn
dfjdr(1) = fe*dftndr/r
!
re    = oxygn%r_e (id_m); 
alpha = oxygn%alph(id_m); A     = oxygn%A   (id_m)
beta  = oxygn%beta(id_m); B     = oxygn%B   (id_m)
kappa = oxygn%kapp(id_m); lambda= oxygn%lamb(id_m)
fe = oxygn%fe; gamm = oxygn%gamm; nu = oxygn%nu

r_re = r/re; r_re_1 = r_re - 1.0D0
lambda20 = (r_re - lambda)**20
lambda19 = lambda20/(r_re - lambda)
kappa20  = (r_re - kappa)**20
kappa19  = kappa20/(r_re - kappa)
exp_alpha = DEXP(-alpha*r_re_1)
exp_beta  = DEXP( -beta*r_re_1)
dphidr = A*exp_alpha*(-alpha*(1.D0 + kappa20) - 20.D0*kappa19) &
     &    /re/(1.D0 + kappa20)**2 - &
     &   B*exp_beta* (-beta*(1.D0 + lambda20) - 20.D0*lambda19) &
     &    /re/(1.D0 + lambda20)**2
!
! ************************* Pair energy and force
phi =  A*exp_alpha/(1.D0 + kappa20) - B*exp_beta/(1.D0 + lambda20)
df  = - dphidr/r
!
! ********************** oxygen part (electron density only)
exp_g  = DEXP(-gamm*r_re_1)
nu20   = (r_re - nu)**20; nu19 = nu20 /(r_re - nu)
fj(2) = fe*exp_g/(1.D0 + nu20)
dfjdr(2) = fe*exp_g*(-gamm*(1.D0 + nu20) - 20.D0*nu19)/re/(1.D0 + nu20)**2/r
!
RETURN
END SUBROUTINE METAL_OXYGEN
!
!#############################################################################
SUBROUTINE OXYGEN_METAL(oxygn, eam, id_m, r, fj, dfjdr, phi, df)
USE DATAFMT
IMPLICIT NONE
TYPE(METL) :: eam
TYPE(OXID) :: oxygn
REAL  (DP) :: r, fj(2), dfjdr(2), phi, df
INTEGER(GI):: id_m
!
REAL   (DP):: re, fe, alpha, beta, lambda, kappa, A, B
REAL   (DP):: r_re, r_re_1, ftn, dftndr, dphidr, gamm, nu, nu20, nu19
REAL   (DP):: lambda20, lambda19, kappa20, kappa19, exp_alpha, exp_beta, exp_g
!
re     = eam%r_e;   fe    = eam%f_e;   beta  = eam%beta; 
lambda = eam%lambda;B     = eam%B
!
r_re = r/re; r_re_1 = r_re - 1.D0
lambda20 = (r_re - lambda)**20; lambda19 = lambda20/(r_re - lambda)
exp_beta  = DEXP(-beta*r_re_1)
ftn  = exp_beta / (1.D0 + lambda20)
dftndr = exp_beta*(-beta*(1.D0 + lambda20) - 20.D0*lambda19) &
     &   /re/(1.D0 + lambda20)**2
!
! ********************** METAL part (electron density only)
fj(2) = fe*ftn
dfjdr(2) = fe*dftndr/r
!
re    = oxygn%r_e (id_m); 
alpha = oxygn%alph(id_m); A     = oxygn%A   (id_m)
beta  = oxygn%beta(id_m); B     = oxygn%B   (id_m)
kappa = oxygn%kapp(id_m); lambda= oxygn%lamb(id_m)
fe = oxygn%fe; gamm = oxygn%gamm; nu = oxygn%nu

r_re = r/re; r_re_1 = r_re - 1.0D0
lambda20 = (r_re - lambda)**20
lambda19 = lambda20/(r_re - lambda)
kappa20  = (r_re - kappa)**20
kappa19  = kappa20/(r_re - kappa)
exp_alpha = DEXP(-alpha*r_re_1)
exp_beta  = DEXP( -beta*r_re_1)
dphidr = A*exp_alpha*(-alpha*(1.D0 + kappa20) - 20.D0*kappa19) &
     &    /re/(1.D0 + kappa20)**2 - &
     &   B*exp_beta* (-beta*(1.D0 + lambda20) - 20.D0*lambda19) &
     &    /re/(1.D0 + lambda20)**2
!
! *********************************** Pair energy and force
phi =  A*exp_alpha/(1.D0 + kappa20) - B*exp_beta/(1.D0 + lambda20)
df  = - dphidr/r
!
! *********************************** Oxygen part (electron density only)
exp_g  = DEXP(-gamm*r_re_1)
nu20   = (r_re - nu)**20; nu19 = nu20 /(r_re - nu)
fj(1) = fe*exp_g/(1.D0 + nu20)
dfjdr(1) = fe*exp_g*(-gamm*(1.D0 + nu20) - 20.D0*nu19)/re/(1.D0 + nu20)**2/r
!
RETURN
END SUBROUTINE OXYGEN_METAL
!
!#############################################################################
SUBROUTINE OXYGEN_ONLY(oxygn, r, fj, dfjdr, phi, df)
USE DATAFMT
IMPLICIT NONE
TYPE(OXID) :: oxygn
REAL  (DP) :: r, fj(2), dfjdr(2), phi, df
!
REAL   (DP):: re, fe, alpha, beta, lambda, kappa, A, B
REAL   (DP):: r_re, r_re_1, dphidr, gamm, nu, nu20, nu19
REAL   (DP):: lambda20, lambda19, kappa20, kappa19, exp_alpha, exp_beta, exp_g
!
!
re    = oxygn%r_e (4); 
alpha = oxygn%alph(4); A     = oxygn%A   (4)
beta  = oxygn%beta(4); B     = oxygn%B   (4)
kappa = oxygn%kapp(4); lambda= oxygn%lamb(4)
fe = oxygn%fe; gamm = oxygn%gamm; nu = oxygn%nu

r_re = r/re; r_re_1 = r_re - 1.0D0
lambda20 = (r_re - lambda)**20
lambda19 = lambda20/(r_re - lambda)
kappa20  = (r_re - kappa)**20
kappa19  = kappa20/(r_re - kappa)
exp_alpha = DEXP(-alpha*r_re_1)
exp_beta  = DEXP( -beta*r_re_1)
dphidr = A*exp_alpha*(-alpha*(1.D0 + kappa20) - 20.D0*kappa19) &
     &    /re/(1.D0 + kappa20)**2 - &
     &   B*exp_beta* (-beta*(1.D0 + lambda20) - 20.D0*lambda19) &
     &    /re/(1.D0 + lambda20)**2
!
! *************************************************** Pair energy and force
phi =  A*exp_alpha/(1.D0 + kappa20) - B*exp_beta/(1.D0 + lambda20)
df  = - dphidr/r
!
! *************************************************** Electron density
exp_g  = DEXP(-gamm*r_re_1)
nu20   = (r_re - nu)**20; nu19 = nu20 /(r_re - nu)
fj(1) = fe*exp_g/(1.D0 + nu20)
fj(2) = fj(1)
dfjdr(1) = fe*exp_g*(-gamm*(1.D0 + nu20) - 20.D0*nu19)/re/(1.D0 + nu20)**2/r
dfjdr(2) = dfjdr(1)
!
RETURN
END SUBROUTINE OXYGEN_ONLY
