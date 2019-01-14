













! DO NOT EDIT!!!

! This file is automatically created by set_variables.py.  To update
! or add variable indices, please edit _variables and then rerun the
! script.


subroutine check_equal(index1, index2)

  use amrex_error_module

  implicit none

  integer, intent(in) :: index1, index2

  if (index1 /= index2) then
    call amrex_error("ERROR: mismatch of indices")
  endif

end subroutine check_equal


subroutine ca_set_auxiliary_indices()


  use meth_params_module
  use network, only: naux, nspec
  implicit none

  QGAMC = 1

  QC = 2

  QDPDR = 3

  QDPDE = 4

end subroutine ca_set_auxiliary_indices

subroutine ca_set_conserved_indices( &
                                    Density, &
                                    Xmom, &
                                    Ymom, &
                                    Zmom, &
                                    Eden, &
                                    Eint, &
                                    Temp, &
                                    FirstAdv, &
                                    FirstSpec, &
                                    FirstAux &
                                   )


  use meth_params_module
  use network, only: naux, nspec
  implicit none
  integer, intent(in) :: Density
  integer, intent(in) :: Xmom
  integer, intent(in) :: Ymom
  integer, intent(in) :: Zmom
  integer, intent(in) :: Eden
  integer, intent(in) :: Eint
  integer, intent(in) :: Temp
  integer, intent(in) :: FirstAdv
  integer, intent(in) :: FirstSpec
  integer, intent(in) :: FirstAux

  URHO = 1
  call check_equal(URHO,Density+1)

  UMX = 2
  call check_equal(UMX,Xmom+1)

  UMY = 3
  call check_equal(UMY,Ymom+1)

  UMZ = 4
  call check_equal(UMZ,Zmom+1)

  UEDEN = 5
  call check_equal(UEDEN,Eden+1)

  UEINT = 6
  call check_equal(UEINT,Eint+1)

  UTEMP = 7
  call check_equal(UTEMP,Temp+1)

  if (nadv > 0) then
    UFA = 8
  else
    UFA = 0
  endif
  call check_equal(UFA,FirstAdv+1)

  if (nspec > 0) then
    UFS = 8 + nadv
  else
    UFS = 0
  endif
  call check_equal(UFS,FirstSpec+1)

  if (naux > 0) then
    UFX = 8 + nadv + nspec
  else
    UFX = 0
  endif
  call check_equal(UFX,FirstAux+1)

end subroutine ca_set_conserved_indices

subroutine ca_set_godunov_indices()


  use meth_params_module
  use network, only: naux, nspec
  implicit none

  GDRHO = 1

  GDU = 2

  GDV = 3

  GDW = 4

  GDPRES = 5

  GDGAME = 6

end subroutine ca_set_godunov_indices

subroutine ca_set_primitive_indices()


  use meth_params_module
  use network, only: naux, nspec
  implicit none

  QRHO = 1

  QU = 2

  QV = 3

  QW = 4

  QGAME = 5

  QPRES = 6

  QREINT = 7

  QTEMP = 8

  QFA = 9

  QFS = 9 + nadv

  QFX = 9 + nadv + nspec

end subroutine ca_set_primitive_indices

