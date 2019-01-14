












module jacobian_sparsity_module


  use burn_type_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains



  subroutine set_jac_entry(state, row, col, val)

    !$acc routine seq

    implicit none

    type (burn_t), intent(inout) :: state
    integer, intent(in) :: row, col
    real(rt), intent(in) :: val

    !$gpu

    state % jac(row, col) = val

  end subroutine set_jac_entry


  subroutine scale_jac_entry(state, row, col, val)

    !$acc routine seq

    implicit none

    type (burn_t), intent(inout) :: state
    integer, intent(in) :: row, col
    real(rt), intent(in) :: val

    !$gpu

    state % jac(row, col) = state % jac(row, col) * val

  end subroutine scale_jac_entry  


  subroutine get_jac_entry(state, row, col, val)

    !$acc routine seq

    implicit none

    type (burn_t), intent(in) :: state
    integer, intent(in) :: row, col
    real(rt), intent(out) :: val

    integer :: csr_loc

    !$gpu

    val = state % jac(row, col)

  end subroutine get_jac_entry


  subroutine set_jac_zero(state)

    !$acc routine seq

    use amrex_constants_module, only: ZERO

    implicit none

    type (burn_t), intent(inout) :: state

    !$gpu

    state % jac(:,:) = ZERO

  end subroutine set_jac_zero

end module jacobian_sparsity_module

