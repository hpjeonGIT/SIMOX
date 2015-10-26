PROGRAM SIMOX !****************************************************************
!
!*** Note
! The code is using the feature of F95/2003. Some of F90 features may not work.
! Required external library
! - FFT for particle mesh Ewald
! - LAPACK for DGESV, linear solver
!
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT) :: Nset
TYPE(TCKT) :: tag
TYPE(METL) :: eam(Nmetal)
TYPE(POTN) :: param(Natom)
TYPE(OXID) :: oxygn
TYPE(STAT) :: sys
TYPE(CORR) :: pair
REAL(DP)   :: dt
INTEGER(GI):: Nloop
TYPE(PTCL), POINTER:: q(:)
TYPE(SORT), POINTER:: cell(:)
REAL(DP),   POINTER:: psi(:,:,:)
REAL(SP)   :: secnds, time0, time1, time2
time0 = 0
time1 = secnds(time0)
!
CALL RANDOM_SEED
CALL READ_INPUT(Nset, tag, sys, q, cell, psi, pair, dt)
CALL CELL_HIERARCHY(Nset, cell)
CALL CELL_SORT(Nset, cell, q)
CALL POTENTIAL_SETUP(eam, param, oxygn)
CALL PME_SELF(Nset, sys, param, q, psi)
IF (.NOT.tag%restart) CALL VINIT(Nset, eam, param, oxygn, sys, q, cell, psi, dt)
Nloop = 1
OPEN(UNIT=15, FILE="temporal_status.dat")
WRITE(15,*) "# time(fs) - T(K) - E_total -   E_Ewald   -   E_eam  -   ", &
     & "Epair   -  Egroup (eV)"
!
IF (tag%regular_run) THEN
   DO WHILE (Nloop <= Nset%freq_max)
      CALL TAG_SET(Nset, tag, Nloop)      
      IF (tag%thermostat .AND. tag%thermo_stoch) THEN 
         CALL VV_NVT1_STOCH(Nset, param, q, dt)
      ELSE; CALL VV_NVE1(Nset, param, q, dt)
      END IF
      IF (tag%cell_sort)  CALL CELL_SORT(Nset, cell, q)
      IF (tag%new_charge) CALL CTIP_SOLVER(Nset, param, sys, q, cell, psi)
      CALL FORCE(Nset, eam, oxygn, sys, param, pair, q, cell)
      CALL PME(Nset, sys, q, psi)
      IF (tag%thermostat) THEN; 
         IF (tag%thermo_stoch) THEN; CALL VV_NVT2_STOCH(Nset, param, sys, q, dt)
         ELSE;  CALL VV_NVT2_BEREND(Nset,param,sys,q,dt);
         END IF
      ELSE; CALL VV_NVE2(Nset,param, sys, q, dt); 
      END IF
      IF (tag%snapshot) CALL PRINT_SNAPSHOT(Nset, Nloop, q)   
      IF (tag%status) CALL PRINT_STAT(Nset, sys, Nloop, dt)
      IF (tag%print_corr) CALL PRINT_RDF(Nset, Nloop, pair)
      Nloop = Nloop + 1
   END DO
ELSE
   DO WHILE (Nloop <= Nset%freq_max)
      CALL TAG_SET(Nset, tag, Nloop)
      CALL VV_NVE1_REFLECT(Nset, param, q, dt); 
      IF (tag%cell_sort) CALL CELL_SORT(Nset, cell, q)
      IF (tag%new_charge) THEN
         CALL REMOVE_RIGID_MOTION(Nset, q)
         CALL OXIDE_CHECK(Nset, tag, q)
         IF (tag%new_o2) THEN
            CALL NEW_O2(Nset, sys, q)
         END IF
         CALL CTIP_SOLVER(Nset, param, sys, q, cell, psi)
      END IF
      CALL FORCE(Nset, eam, oxygn, sys, param, pair, q, cell)
      CALL PME(Nset, sys, q, psi)
      IF (tag%thermostat) THEN; CALL VV_NVT2_BEREND(Nset, param, sys, q, dt); 
      ELSE; CALL VV_NVE2(Nset,param, sys, q, dt); 
      END IF
      IF (tag%snapshot) CALL PRINT_SNAPSHOT(Nset, Nloop, q)
      IF (tag%status) CALL PRINT_STAT(Nset, sys, Nloop, dt)
      IF (tag%print_corr) CALL PRINT_RDF(Nset, Nloop, pair)
      Nloop = Nloop + 1
   END DO
END IF
!
CALL EVACUATE_MEMORY(q, psi)
CLOSE(15)
time2 = secnds(time1)
PRINT '("# Wall time =", ES11.4, " seconds with", I8, " loops")', time2, Nloop-1
STOP
!
CONTAINS
!
!###############################################################################
SUBROUTINE READ_INPUT(Nset, tag, sys, q, cell, psi, pair, dt)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(TCKT):: tag
TYPE(STAT):: sys
TYPE(CORR):: pair
REAL(DP)  :: dt
TYPE(PTCL), POINTER:: q(:)
TYPE(SORT), POINTER:: cell(:)
REAL(DP),   POINTER:: psi(:,:,:)
!
INTEGER(GI):: openstatus, istatus, i, id, Ncell(3)
REAL   (DP):: Lcell(3)
CHARACTER(256):: dummy, FILENAME
!
! Initialize tag variables
tag%regular_run = .FALSE.
tag%new_o2      = .FALSE.
tag%new_charge  = .FALSE.
tag%snapshot    = .FALSE.
tag%thermostat  = .FALSE.
tag%restart     = .FALSE.
tag%status      = .FALSE.
!
! Read simulation parameters **************************************************
OPEN(UNIT=10, FILE='config.prm', STATUS = "OLD", ACTION="READ", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) STOP "#### cannot  open config.prm ###"
READ(10,*) dummy
READ(10,*) dummy
READ(10,*) dummy
CALL ANY_2_UPPER(dummy)
IF (dummy == 'REGULAR') THEN
   tag%regular_run = .TRUE.
ELSE IF (dummy == 'OXIDATION') THEN
   tag%regular_run = .FALSE.
ELSE
   STOP "=== Simulation type is unknown ==="
END IF
READ(10,*) dummy
READ(10,*) sys%time_max, sys%time_snap, sys%time_rdf, sys%time_stat,dt
Nset%freq_max  = NINT(sys%time_max /dt)
Nset%freq_snap = NINT(sys%time_snap/dt)
Nset%freq_rdf  = NINT(sys%time_rdf /dt)
Nset%freq_stat = NINT(sys%time_stat /dt)
dt = dt/ TFM
READ(10,*) dummy
READ(10,*) sys%box(:)
Nset%box(:) = sys%box(:)
Nset%V = Nset%box(1)*Nset%box(2)*Nset%box(3)
READ(10,*) dummy
READ(10,*) Nset%fft(:)
ALLOCATE(psi(Nset%fft(1),Nset%fft(2),Nset%fft(3)), STAT=istatus)
READ(10,*) dummy
READ(10,*) dummy, sys%T_given, Nset%alpha
sys%T_given = sys%T_given*kB ! K -> eV
CALL ANY_2_UPPER(dummy)
IF (dummy == 'YES') THEN
   tag%thermostat = .TRUE.
ELSE IF (dummy == 'NO') THEN
   tag%thermostat = .FALSE.
ELSE
   STOP "=== Temperature control is unknown ==="
END IF
READ(10,*) dummy
READ(10,*) dummy
CALL ANY_2_UPPER(dummy)
IF (dummy == 'YES') THEN
   tag%restart = .TRUE.
ELSE IF (dummy == 'NO') THEN
   tag%restart = .FALSE.
ELSE
   STOP "=== Restart option is unknown ==="
END IF
READ(10,*) dummy
READ(10,*) sys%Rcut
sys%Rcut2 = sys%Rcut**2
sys%a = 3.8D0/sys%Rcut
READ(10,*) dummy
READ(10,*) Nset%freq_new, Nset%freq_sort
READ(10,*) dummy
READ(10,*) Nset%z_new, Nset%N_o2
READ(10,*) dummy
READ(10,*) Nset%q_crit
CLOSE(10)
IF (Nset%freq_new < 0) THEN
   tag%thermo_stoch = .TRUE.
ELSE
   tag%thermo_stoch = .FALSE.
END IF
!
! Define cell (linked list set) ***********************************************
Ncell(:) = INT(Nset%box(:)/sys%Rcut)
Lcell(:) = Nset%box(:)/DBLE(Ncell(:))
IF (Ncell(1) < 3 .AND. Ncell(2) < 3 .AND. Ncell(3) < 3) THEN
   PRINT *, "==== Error - unitcell is too small or cut-off radius is too small"
   STOP
END IF
WRITE(*,50) Ncell(:), Lcell(:)
50 FORMAT("# === ",I3,"x",I3,"x",I3," cells with size of", 3(F5.1,1X), "===")
Nset%Ncell(:)  = Ncell(:)
Nset%Ncell_all = Ncell(1)*Ncell(2)*Ncell(3)
Nset%Lcell(:)  = Lcell(:)
!
ALLOCATE(cell(Nset%Ncell_all), STAT=istatus)
IF(istatus /=0) STOP "=== Cell allocation error ==="
!
! Read particle data **********************************************************
IF (tag%restart) THEN
   FILENAME = 'restart.xyz'
ELSE
   FILENAME = 'input.xyz'
END IF
OPEN(UNIT=11, FILE=FILENAME, STATUS = "OLD", ACTION="READ", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) THEN
   PRINT *, '("=== cannot open ", A, "===")', FILENAME; STOP
END IF
!
! Atomic data configuration
! Al=1 Ni=2 Fe=3 O=4 H=5
Nset%Npt(:) = 0
READ(11,*) Nset%Npt_all
ALLOCATE(q(Nset%Npt_all), STAT=istatus)
IF(istatus /=0) STOP "=== Particle allocation error ==="
READ(11,*) dummy
DO i=1, Nset%Npt_all
   IF (tag%restart) THEN
      READ(11,*) dummy, q(i)%xx(:), q(i)%xv(:), q(i)%ff(:), q(i)%q
   ELSE
      READ(11,*) dummy, q(i)%xx(:), q(i)%q
   END IF
   SELECT CASE (dummy)
   CASE ("Al")
      id = 1
   CASE ("Ni")
      id = 2
   CASE ("Fe")
      id = 3
   CASE ("O")
      id = 4
   CASE ("H")
      id = 5
   CASE DEFAULT
      STOP "=== Particle type error ==="
   END SELECT
   Nset%Npt(id) = Nset%Npt(id) + 1
   q(i)%id = id
END DO
CLOSE(11)
id = 0; DO i=1,Natom; id = id + Nset%Npt(i); END DO
IF (id /= Nset%Npt_all) THEN
   DEALLOCATE(q, STAT = istatus)
   STOP "=== Number of particles mismatch ==="
END IF
!
! Correlation data setup ******************************************************
pair%dr = sys%Rcut / DBLE(Nrdf)
pair%V  = Nset%V
pair%rdf(:,:,:) = 0
pair%Npt_old    = Nset%Npt_all
!
RETURN 
END SUBROUTINE READ_INPUT
!
!##############################################################################
SUBROUTINE EVACUATE_MEMORY(q, psi)
USE DATAFMT
IMPLICIT NONE
!
TYPE(PTCL), POINTER:: q(:)
REAL(DP),   POINTER:: psi(:,:,:)
INTEGER(GI):: istatus
DEALLOCATE(q, STAT = istatus)
IF(istatus /=0) STOP "=== Particle deallocation error ==="
DEALLOCATE(psi, STAT = istatus)
IF(istatus /=0) STOP "=== PME coefficient deallocation error ==="
!
RETURN
END SUBROUTINE EVACUATE_MEMORY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END PROGRAM SIMOX
!
!http://coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2005-03/0761.html
!###############################################################################
SUBROUTINE ANY_2_UPPER(txt_string)
USE DATAFMT, ONLY: GI
IMPLICIT NONE
CHARACTER(LEN=*):: txt_string
INTEGER(GI):: i, nlen, id
nlen = LEN(txt_string)
DO i=1, nlen
   id = ichar(txt_string(i:i))
   IF (id >= 97 .AND. id < 122) txt_string(i:i) = CHAR(id-32)
END DO
RETURN
END SUBROUTINE ANY_2_UPPER
!
!###############################################################################
SUBROUTINE TAG_SET(Nset, tag, Nloop)
USE DATAFMT, ONLY: GI, DP, AMNT, TCKT, STAT
IMPLICIT NONE
!
TYPE(AMNT) :: Nset
TYPE(TCKT) :: tag
INTEGER(GI):: Nloop
!
IF (Nset%freq_new < 0) THEN
   tag%new_charge = .FALSE.
ELSE
   IF (MOD(Nloop, Nset%freq_new) == 1) THEN
      tag%new_charge = .TRUE.
   ELSE 
      tag%new_charge = .FALSE.
   END IF
END IF
IF (MOD(Nloop, Nset%freq_snap) == 0) THEN
   tag%snapshot = .TRUE.
ELSE
   tag%snapshot = .FALSE.
END IF
IF (MOD(Nloop, Nset%freq_stat) == 0) THEN
   tag%status = .TRUE.
ELSE
   tag%status = .FALSE.
END IF
IF (MOD(Nloop, Nset%freq_rdf) == 0) THEN
   tag%print_corr = .TRUE.
ELSE
   tag%print_corr = .FALSE.
END IF
IF (MOD(Nloop, Nset%freq_sort) == 0) THEN
   tag%cell_sort = .TRUE.
ELSE
   tag%cell_sort = .FALSE.
END IF
RETURN
END SUBROUTINE TAG_SET
