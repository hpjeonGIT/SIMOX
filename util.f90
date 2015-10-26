!
!#############################################################################
SUBROUTINE PRINT_STAT(Nset, sys, Nloop, dt)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(STAT):: sys
INTEGER(GI):: Nloop
REAL(DP):: dt
!
REAL(DP):: E_total, E_es, E_eam, Temperature
E_es = sys%Eself + sys%Edir + sys%Erec + sys%Ectip
E_eam = sys%Epair + sys%Egroup
E_total = E_es + E_eam + sys%mv2/2.D0
Temperature = sys%mv2/DBLE(Nset%Npt_all*3)/kB
WRITE(15,100) Nloop*TFM*dt, Temperature, E_total, E_es, E_eam, sys%Epair, &
     & sys%Egroup
100 FORMAT (ES11.4, 1X, F6.1, 5(1X, ES11.4))
RETURN
END SUBROUTINE PRINT_STAT
!
!#############################################################################
SUBROUTINE PRINT_SNAPSHOT(Nset, Nloop, q)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(PTCL):: q(Nset%Npt_all)
INTEGER(GI):: Nloop
!
INTEGER:: openstatus, nframe, i, id
CHARACTER(LEN=256):: dummy, RESTFILE, SNAPFILE
!
nframe = Nloop/Nset%freq_snap
WRITE(RESTFILE,100) nframe
100 FORMAT("snapshot", I3.3, ".xyz")
WRITE(SNAPFILE,110) nframe
110 FORMAT("image", I3.3, ".xyz")
!
OPEN(UNIT=20, FILE=RESTFILE, ACTION="WRITE", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) THEN
   PRINT '("=== cannot open ", A, "===")', RESTFILE ; STOP
END IF
OPEN(UNIT=21, FILE=SNAPFILE, ACTION="WRITE", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) THEN
   PRINT '("=== cannot open ", A, "===")', SNAPFILE ; STOP
END IF

WRITE(20,*) Nset%Npt_all
WRITE(20,*) "Frame = ", Nframe, " energy =  0.0 "
WRITE(21,*) Nset%Npt_all
WRITE(21,*) "Frame = ", Nframe, " energy =  0.0 "

DO i=1, Nset%Npt_all
   id = q(i)%id
   SELECT CASE (id)
      CASE (1)
         dummy = "Al"
      CASE (2)
         dummy = "Ni"
      CASE (3)
         dummy = "Fe"
      CASE(4)
         dummy = "O "
      CASE(5)
         dummy = "H "
      CASE DEFAULT
         CLOSE(20)
         PRINT *, id
         STOP "=== Particle id error ==="
   END SELECT
   WRITE(20,200) dummy, q(i)%xx(:), q(i)%xv(:), q(i)%ff(:), q(i)%q
   WRITE(21,300) dummy, q(i)%xx(:), q(i)%q
END DO
200 FORMAT(A3, 3(ES13.6, 1X), 6(ES11.4, 1X), ES13.6)
300 FORMAT(A3, 4(ES11.4, 1X))
CLOSE(20)
CLOSE(21)
RETURN
END SUBROUTINE PRINT_SNAPSHOT
!
!****************************************************************************** 
FUNCTION fluct(x)
  USE DATAFMT
  IMPLICIT NONE
  REAL(DP):: fluct, x, r, v1, v2
  REAL(DP):: rand1, rand2
  !
  ! Initialization
  r=1.D0
  DO WHILE (r.ge.1.D0)
     CALL RANDOM_NUMBER(rand1)
     CALL RANDOM_NUMBER(rand2)
     v1 = 2.D0*rand1 - 1.D0
     v2 = 2.D0*rand2 - 1.D0
     r = v1*v1+v2*v2
  END DO
  fluct = v1*DSQRT(-2.D0*DLOG(r)/r)*x
  RETURN
END FUNCTION fluct
!
!###############################################################################
SUBROUTINE RANDOM_POSITION_O2(n_o, z_new, Lx, Ly, xx)
USE DATAFMT, ONLY: DP, GI, PI
INTEGER(GI):: n_o
REAL   (DP):: z_new, Lx, Ly, xx(n_o, 3)
!
INTEGER(GI):: i, j, n_o2
REAL   (DP):: cx(n_o/2, 3), rand, phi, theta, L0, Rcut, Rcut2, r2
LOGICAL    :: pass
!
L0 = 2.0D0 ! distance betwen O-O
n_o2 = n_o /2 
Rcut = DSQRT(Lx*Ly/DBLE(n_o2))*0.5D0
IF (Rcut < L0) STOP "=== Too many new oxygens ==="
Rcut2 = Rcut**2
!
CALL RANDOM_NUMBER(rand);   cx(1,1) = (rand - 0.5D0)*Lx
CALL RANDOM_NUMBER(rand);   cx(1,2) = (rand - 0.5D0)*Ly
DO i = 2, n_o2
   pass = .TRUE.
   DO WHILE (pass)
      pass = .FALSE.
      CALL RANDOM_NUMBER(rand);   cx(i,1) = (rand - 0.5D0)*Lx
      CALL RANDOM_NUMBER(rand);   cx(i,2) = (rand - 0.5D0)*Ly
      DO j = 1, i-1
         r2 = (cx(i,1) - cx(j,1))**2 + (cx(i,2) - cx(j,2))**2
         IF (r2 < Rcut2) pass = .TRUE.
      END DO
   END DO
END DO
cx(:,3) = z_new
!
DO i = 1, n_o2
   CALL RANDOM_NUMBER(rand); phi = rand*PI*2.D0
   CALL RANDOM_NUMBER(rand); theta = rand*PI
   xx(i*2-1,1) = L0*DSIN(theta)*DCOS(phi)
   xx(i*2-1,2) = L0*DSIN(theta)*DCOS(phi)
   xx(i*2-1,3) = L0*DCOS(theta)
   xx(i*2  ,1) = - xx(i*2-1,1)
   xx(i*2  ,2) = - xx(i*2-1,2)
   xx(i*2  ,3) = - xx(i*2-1,3)
   xx(i*2-1,:) = xx(i*2-1,:) + cx(i,:)
   xx(i*2  ,:) = xx(i*2  ,:) + cx(i,:)
END DO
!
RETURN
END SUBROUTINE RANDOM_POSITION_O2
!
!#############################################################################
SUBROUTINE OXIDE_CHECK(Nset, tag, q)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(TCKT):: tag
TYPE(PTCL):: q(Nset%Npt_all)
!
INTEGER(GI):: Npt_all, i
REAL   (DP):: o_height, m_height, q_max
!
! Highest metal position
Npt_all = Nset%Npt_all
m_height = 0.0D0
DO i=1, Npt_all
   IF (q(i)%id == 1) THEN ! Al. only
      m_height = MAX(q(i)%xx(3), m_height)
   END IF
END DO
!
! Only newly inserted O2
o_height = 0.0D0
q_max = -10.0D0
DO i=Npt_all, Npt_all-Nset%N_o2*2 + 1, -1
   IF (q(i)%id == 4) THEN
      o_height = MAX(q(i)%xx(3), o_height)
      q_max    = MAX(q(i)%q, q_max)
   END IF
END DO
!
IF ((m_height + 0.0D0) > o_height .OR. q_max < Nset%q_crit) THEN
!IF ( q_max < Nset%q_crit) THEN
   tag%new_o2 = .TRUE.
ELSE
   tag%new_o2 = .FALSE.
END IF
RETURN
END SUBROUTINE OXIDE_CHECK
!
!#############################################################################
SUBROUTINE OXIDE_CHECK2(Nset, tag, q)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(TCKT):: tag
TYPE(PTCL):: q(Nset%Npt_all)
!
INTEGER(GI):: Npt_all, i, Npt_test
!
Npt_all = Nset%Npt_all
Npt_test = Npt_all - Nset%N_o2*2 + 1
tag%new_o2 = .TRUE.
DO i=Npt_all, Npt_test, -1
   IF (q(i)%id == 4 .AND. q(i)%q > Nset%q_crit) tag%new_o2 = .FALSE.
END DO
!
END SUBROUTINE OXIDE_CHECK2
!
!#############################################################################
SUBROUTINE OXIDE_CHECK_OLD(Nset, tag, q)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
TYPE(TCKT):: tag
TYPE(PTCL):: q(Nset%Npt_all)
!
INTEGER(GI):: Npt_all, i
REAL   (DP):: o_height, m_height
!
! Highest metal position
Npt_all = Nset%Npt_all
m_height = 0.0D0
DO i=1, Npt_all
   IF (q(i)%id < 4) THEN
      m_height = MAX(q(i)%xx(3), m_height)
   END IF
END DO
!
! Only newly inserted O2
o_height = 0.0D0
DO i=Npt_all, Npt_all-Nset%N_o2*2 + 1, -1
   IF (q(i)%id == 4) THEN
      o_height = MAX(q(i)%xx(3), o_height)
   END IF
END DO
!
IF (m_height > o_height) THEN
   tag%new_o2 = .TRUE.
ELSE
   tag%new_o2 = .FALSE.
END IF
RETURN
END SUBROUTINE OXIDE_CHECK_OLD
!
!#############################################################################
SUBROUTINE PRINT_RDF(Nset, Nloop, pair)
USE DATAFMT
IMPLICIT NONE
TYPE(AMNT):: Nset
INTEGER(GI):: Nloop
TYPE(CORR):: pair
!
INTEGER(GI):: openstatus, nframe, i, j, n, Nkind, Npt_old, Npt_new
REAL(DP)   :: result(Natom*(Natom+1)/2, Nrdf), r, dr, const, loop_sum, Xrdf
CHARACTER(LEN=256):: FILENAME

!
nframe = Nloop/Nset%freq_rdf
WRITE(FILENAME,100) nframe
100 FORMAT("rdf", I3.3, ".dat")
!
OPEN(UNIT=20, FILE=FILENAME, ACTION="WRITE", &
     & FORM = "FORMATTED", IOSTAT = openstatus)
IF (openstatus >0) THEN
   PRINT '("=== cannot open ", A, "===")', FILENAME ; STOP
END IF
!
Xrdf = DBLE(Nset%freq_rdf)
n = 0
DO i=1, Natom
   DO j=i, Natom
      n = n + 1
      IF (i /= j) THEN
         result(n, 1:Nrdf) = DBLE(pair%rdf(i,j,1:Nrdf) + pair%rdf(j,i,1:Nrdf)) &
              &  /Xrdf
      ELSE
         result(n, 1:Nrdf) = DBLE(pair%rdf(i,j,1:Nrdf))/Xrdf
      END IF
   END DO
END DO
dr = pair%dr
Nkind = Natom*(Natom+1)/2
Npt_old = pair%Npt_old
Npt_new = pair%Npt_new
loop_sum = DBLE(Npt_old*Npt_new)*4.D0*PI/3.0D0/pair%V
WRITE(20,300) Nloop, Npt_old, Npt_new
300 FORMAT("# Header: current loop = ",I7," Npt_old = ", I5, " Npt_new = ", I5)
WRITE(20,400)
400 FORMAT("# Radius(A) / Al-Al(2) / Al-Ni(3) / Al-Fe(4) / Al-O(5) / Al-H(6) " &
         & "/ Ni-Ni(7) / Ni-Fe(8) / Ni-O(9) / Ni-H(10) / Fe-Fe(11) / Fe-O(12)" &
         & " / Fe-H(13) / O-O(14) / O-H(15) / H-H(16) ")
DO i=1, Nrdf
   r = DBLE(i)*dr
   const = ((r+dr)**3 - dr**3)*loop_sum
   WRITE(20,500, ADVANCE="NO") dr*(DBLE(i)-0.5D0)
   DO j=1, Nkind
      WRITE(20,500, ADVANCE="NO") result(j,i)/const
   END DO
   WRITE(20,*)
END DO
500 FORMAT(ES10.3, 1X)
CLOSE(20)
pair%rdf(:,:,:) = 0
pair%Npt_old = pair%Npt_new
!
RETURN
END SUBROUTINE PRINT_RDF

