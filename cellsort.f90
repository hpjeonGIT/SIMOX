!
!############################################################################
SUBROUTINE CELL_HIERARCHY(Nset, cell)
USE DATAFMT
IMPLICIT NONE
!
TYPE(AMNT):: Nset
TYPE(SORT):: cell(Nset%Ncell_all)
!
INTEGER(GI):: i, j, k, Ncell(3), ID_CELL, id
Ncell(:) = Nset%Ncell(:)
!
DO i=1, Ncell(1)
   DO j=1, Ncell(2)
      DO k=1, Ncell(3)
         id = ID_CELL(Ncell, i, j, k)

         cell(id)%pair( 1) = ID_CELL(Ncell, i-1, j-1, k-1)
         cell(id)%pair( 2) = ID_CELL(Ncell, i  , j-1, k-1)
         cell(id)%pair( 3) = ID_CELL(Ncell, i+1, j-1, k-1)
         cell(id)%pair( 4) = ID_CELL(Ncell, i-1, j  , k-1)
         cell(id)%pair( 5) = ID_CELL(Ncell, i  , j  , k-1)
         cell(id)%pair( 6) = ID_CELL(Ncell, i+1, j  , k-1)
         cell(id)%pair( 7) = ID_CELL(Ncell, i-1, j+1, k-1)
         cell(id)%pair( 8) = ID_CELL(Ncell, i  , j+1, k-1)
         cell(id)%pair( 9) = ID_CELL(Ncell, i+1, j+1, k-1)
         cell(id)%pair(10) = ID_CELL(Ncell, i-1, j-1, k  )
         cell(id)%pair(11) = ID_CELL(Ncell, i  , j-1, k  )
         cell(id)%pair(12) = ID_CELL(Ncell, i+1, j-1, k  )
         cell(id)%pair(13) = ID_CELL(Ncell, i-1, j  , k  )
      END DO
   END DO
END DO
!
RETURN
END SUBROUTINE CELL_HIERARCHY
!
!*****************************************************************************
FUNCTION ID_CELL(Ncell, i, j, k)
USE DATAFMT, ONLY:GI
IMPLICIT NONE
!
INTEGER(GI):: Ncell(3), i, j, k, ID_CELL
INTEGER(GI):: l, m, n
l = MOD(i-1+Ncell(1), Ncell(1))
m = MOD(j-1+Ncell(2), Ncell(2))
n = MOD(k-1+Ncell(3), Ncell(3))
ID_CELL = l + m*Ncell(1) + n*Ncell(1)*Ncell(2) + 1
RETURN
END FUNCTION ID_CELL
!
!############################################################################
SUBROUTINE CELL_SORT(Nset, cell, q)
USE DATAFMT
TYPE(AMNT):: Nset
TYPE(SORT):: cell(Nset%Ncell_all)
TYPE(PTCL):: q(Nset%Npt_all)
!
INTEGER(GI):: i, Npt(Nset%Ncell_all), Npt_all, Ncell(3), id(3), id_cell
REAL   (DP):: hbox(3), Lcell(3), xx(3)
!
Npt_all = Nset%Npt_all
hbox(:) = 0.5D0*Nset%box(:)
Ncell(:) = Nset%Ncell(:)
Lcell(:) = Nset%Lcell(:)
!
Npt(:) = 0
DO i=1, Npt_all
   xx(:) = q(i)%xx(:) + hbox(:)
   id(:) = INT(xx(:)/Lcell(:))
   id(:) = MOD(id(:) + Ncell(:), Ncell(:))
   id_cell = id(1) + id(2)*Ncell(1) + id(3)*Ncell(1)*Ncell(2) + 1
   Npt(id_cell) = Npt(id_cell) + 1
   cell(id_cell)%link(Npt(id_cell)) = i
   IF (Npt(id_cell) > Npt_max) STOP "=== Cell particle limit crashed !!! ==="
END DO
cell(:)%Npt = Npt(:)
!
RETURN
END SUBROUTINE CELL_SORT
