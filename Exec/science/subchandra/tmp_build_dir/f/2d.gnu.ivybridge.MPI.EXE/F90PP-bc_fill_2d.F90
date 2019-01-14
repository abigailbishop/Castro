












module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine hypfill(adv, adv_l1, adv_l2, adv_h1, adv_h2, &
                     domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_h1, adv_h2
    integer,  intent(in   ) :: bc(2,2,NVAR)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    integer  :: lo(3), hi(3)

    lo(1) = adv_l1
    lo(2) = adv_l2
    lo(3) = 0
    hi(1) = adv_h1
    hi(2) = adv_h2
    hi(3) = 0

    call amrex_filccn(lo, hi, adv, lo, hi, NVAR, domlo, domhi, delta, xlo, bc)

  end subroutine hypfill


  subroutine ca_hypfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_hypfill")

    use bc_ext_fill_module, only: ext_fill

    implicit none

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_h1, adv_h2
    integer,  intent(in   ) :: bc(2,2,NVAR)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2,NVAR)

    call hypfill(adv, adv_l1, adv_l2, adv_h1, adv_h2, domlo, domhi, delta, xlo, time, bc)

    ! process the external BCs here
    call ext_fill(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,time,bc)

  end subroutine ca_hypfill


  subroutine denfill(adv, adv_l1, adv_l2, adv_h1, adv_h2, &
                     domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_h1, adv_h2
    integer,  intent(in   ) :: bc(2,2,1)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2)

    integer :: lo(3), hi(3)

    lo(1) = adv_l1
    lo(2) = adv_l2
    lo(3) = 0
    hi(1) = adv_h1
    hi(2) = adv_h2
    hi(3) = 0

    call amrex_filccn(lo, hi, adv, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine denfill


  subroutine ca_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2, &
                        domlo,domhi,delta,xlo,time,bc) bind(C, name="ca_denfill")

    use bc_ext_fill_module, only: ext_denfill

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: adv_l1, adv_l2, adv_h1, adv_h2
    integer,  intent(in   ) :: bc(2,2,1)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: adv(adv_l1:adv_h1,adv_l2:adv_h2)

    call denfill(adv, adv_l1, adv_l2, adv_h1, adv_h2, domlo, domhi, delta, xlo, time, bc)

    ! process the external BCs here
    call ext_denfill(adv,adv_l1,adv_l2,adv_h1,adv_h2,domlo,domhi,delta,xlo,time,bc)

  end subroutine ca_denfill


  subroutine phigravfill(phi, phi_l1, phi_l2, phi_h1, phi_h2, &
                         domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_h1, phi_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    integer :: lo(3), hi(3)

    lo(1) = phi_l1
    lo(2) = phi_l2
    lo(3) = 0
    hi(1) = phi_h1
    hi(2) = phi_h2
    hi(3) = 0

    call amrex_filccn(lo, hi, phi, lo, hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine phigravfill


  subroutine ca_phigravfill(phi,phi_l1,phi_l2,phi_h1,phi_h2, &
                            domlo,domhi,delta,xlo,time,bc) &
                            bind(C, name="ca_phigravfill")

    implicit none

    integer,  intent(in   ) :: phi_l1, phi_l2, phi_h1, phi_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: phi(phi_l1:phi_h1,phi_l2:phi_h2)

    call phigravfill(phi, phi_l1, phi_l2, phi_h1, phi_h2, &
                     domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_phigravfill


  subroutine gravxfill(grav, grav_l1, grav_l2, grav_h1, grav_h2, &
                       domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: lo(3), hi(3)
    integer :: bc_temp(2,2)
    integer :: d

    lo(1) = grav_l1
    lo(2) = grav_l2
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = grav_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here
    bc_temp(:,:) = bc(:,:)

    do d = 1, 2
       if (bc(d,1) == EXT_DIR .and. lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       endif

       if (bc(d,2) == EXT_DIR .and. hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       endif
    end do

    call amrex_filccn(lo, hi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine gravxfill


  subroutine ca_gravxfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravxfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call gravxfill(grav, grav_l1, grav_l2, grav_h1, grav_h2, &
                   domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravxfill


  subroutine gravyfill(grav, grav_l1, grav_l2, grav_h1, grav_h2, &
                       domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: lo(3), hi(3)
    integer :: bc_temp(2,2)
    integer :: d

    lo(1) = grav_l1
    lo(2) = grav_l2
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = grav_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here
    bc_temp(:,:) = bc(:,:)

    do d = 1, 2
       if (bc(d,1) == EXT_DIR .and. lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       endif

       if (bc(d,2) == EXT_DIR .and. hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       endif
    end do

    call amrex_filccn(lo, hi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine gravyfill


  subroutine ca_gravyfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravyfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call gravyfill(grav, grav_l1, grav_l2, grav_h1, grav_h2, &
                   domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravyfill


  subroutine gravzfill(grav, grav_l1, grav_l2, grav_h1, grav_h2, &
                       domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    integer :: lo(3), hi(3)
    integer :: bc_temp(2,2)
    integer :: d

    lo(1) = grav_l1
    lo(2) = grav_l2
    lo(3) = 0
    hi(1) = grav_h1
    hi(2) = grav_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here
    bc_temp(:,:) = bc(:,:)

    do d = 1, 2
       if (bc(d,1) == EXT_DIR .and. lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       endif

       if (bc(d,2) == EXT_DIR .and. hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       endif
    end do

    call amrex_filccn(lo, hi, grav, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine gravzfill


  subroutine ca_gravzfill(grav,grav_l1,grav_l2,grav_h1,grav_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_gravzfill")

    implicit none

    integer,  intent(in   ) :: grav_l1, grav_l2, grav_h1, grav_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: grav(grav_l1:grav_h1,grav_l2:grav_h2)

    call gravzfill(grav, grav_l1, grav_l2, grav_h1, grav_h2, &
                   domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_gravzfill






  subroutine reactfill(react, react_l1, react_l2, react_h1, react_h2, &
                       domlo, domhi, delta, xlo, time, bc)

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: react_l1, react_l2, react_h1, react_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: react(react_l1:react_h1,react_l2:react_h2)

    integer :: lo(3), hi(3)
    integer :: bc_temp(2,2)
    integer :: d

    lo(1) = react_l1
    lo(2) = react_l2
    lo(3) = 0
    hi(1) = react_h1
    hi(2) = react_h2
    hi(3) = 0

    ! handle an external BC via extrapolation here
    bc_temp(:,:) = bc(:,:)

    do d = 1, 2
       if (bc(d,1) == EXT_DIR .and. lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       endif

       if (bc(d,2) == EXT_DIR .and. hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       endif
    end do

    call amrex_filccn(lo, hi, react, lo, hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine reactfill


  subroutine ca_reactfill(react,react_l1,react_l2,react_h1,react_h2, &
                          domlo,domhi,delta,xlo,time,bc) &
                          bind(C, name="ca_reactfill")

    implicit none

    integer,  intent(in   ) :: react_l1, react_l2, react_h1, react_h2
    integer,  intent(in   ) :: bc(2,2)
    integer,  intent(in   ) :: domlo(2), domhi(2)
    real(rt), intent(in   ) :: delta(2), xlo(2), time
    real(rt), intent(inout) :: react(react_l1:react_h1,react_l2:react_h2)

    call reactfill &
         (react, react_l1, react_l2, react_h1, react_h2, &
         domlo, domhi, delta, xlo, time, bc)

  end subroutine ca_reactfill




end module bc_fill_module
