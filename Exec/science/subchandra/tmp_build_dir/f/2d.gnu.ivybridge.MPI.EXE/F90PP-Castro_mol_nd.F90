












! advection routines in support of method of lines integration

subroutine ca_mol_single_stage(lo, hi, time, &
                               domlo, domhi, &
                               stage_weight, &
                               uin, uin_lo, uin_hi, &
                               uout, uout_lo, uout_hi, &
                               q, q_lo, q_hi, &
                               qaux, qa_lo, qa_hi, &
                               srcU, srU_lo, srU_hi, &
                               update, updt_lo, updt_hi, &
                               update_flux, uf_lo, uf_hi, &
                               dx, dt, &
                               flux1, flux1_lo, flux1_hi, &
                               flux2, flux2_lo, flux2_hi, &
                               area1, area1_lo, area1_hi, &
                               area2, area2_lo, area2_hi, &
                               pradial, p_lo, p_hi, &
                               dloga, dloga_lo, dloga_hi, &
                               vol, vol_lo, vol_hi, &
                               verbose) bind(C, name="ca_mol_single_stage")

  use amrex_error_module
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use meth_params_module, only : NQ, QVAR, NVAR, NGDNV, GDPRES, &
                                 UTEMP, UEINT, USHK, GDU, GDV, GDW, UMX, &
                                 use_flattening, QPRES, NQAUX, &
                                 QTEMP, QFS, QFX, QREINT, QRHO, &
                                 first_order_hydro, difmag, hybrid_riemann, &
                                 limit_fluxes_on_small_dens, ppm_type, ppm_temp_fix
  use advection_util_module, only : limit_hydro_fluxes_on_small_dens, shock, &
                                    divu, normalize_species_fluxes, calc_pdivu, &
                                    scale_flux, apply_av
  use amrex_constants_module, only : ZERO, HALF, ONE, FOURTH
  use flatten_module, only: uflatten
  use riemann_module, only: cmpflx
  use ppm_module, only : ppm_reconstruct
  use amrex_fort_module, only : rt => amrex_real
  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module, only : eos
  use network, only : nspec, naux
  use prob_params_module, only : dg, coord_type


  implicit none

  integer, intent(in) :: lo(3), hi(3), verbose
  integer, intent(in) ::  domlo(3), domhi(3)
  real(rt), intent(in) :: stage_weight
  integer, intent(in) :: uin_lo(3), uin_hi(3)
  integer, intent(in) :: uout_lo(3), uout_hi(3)
  integer, intent(in) :: q_lo(3), q_hi(3)
  integer, intent(in) :: qa_lo(3), qa_hi(3)
  integer, intent(in) :: srU_lo(3), srU_hi(3)
  integer, intent(in) :: updt_lo(3), updt_hi(3)
  integer, intent(in) :: uf_lo(3), uf_hi(3)
  integer, intent(in) :: flux1_lo(3), flux1_hi(3)
  integer, intent(in) :: area1_lo(3), area1_hi(3)
  integer, intent(in) :: flux2_lo(3), flux2_hi(3)
  integer, intent(in) :: area2_lo(3), area2_hi(3)
  integer, intent(in) :: p_lo(3), p_hi(3)
  integer, intent(in) :: dloga_lo(3), dloga_hi(3)
  integer, intent(in) :: vol_lo(3), vol_hi(3)

  real(rt), intent(in) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
  real(rt), intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
  real(rt), intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
  real(rt), intent(inout) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
  real(rt), intent(in) :: srcU(srU_lo(1):srU_hi(1), srU_lo(2):srU_hi(2), srU_lo(3):srU_hi(3), NVAR)
  real(rt), intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
  real(rt), intent(inout) :: update_flux(uf_lo(1):uf_hi(1), uf_lo(2):uf_hi(2), uf_lo(3):uf_hi(3), NVAR)
  real(rt), intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
  real(rt), intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
  real(rt), intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
  real(rt), intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
  real(rt), intent(inout) :: pradial(p_lo(1):p_hi(1), p_lo(2):p_hi(2), p_lo(3):p_hi(3))
  real(rt), intent(in) :: dloga(dloga_lo(1):dloga_hi(1), dloga_lo(2):dloga_hi(2), dloga_lo(3):dloga_hi(3))
  real(rt), intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))
  real(rt), intent(in) :: dx(3), dt, time

  ! Automatic arrays for workspace
  real(rt)        , pointer:: flatn(:,:,:)
  real(rt)        , pointer:: div(:,:,:)

  ! Edge-centered primitive variables (Riemann state)
  real(rt)        , pointer:: q1(:,:,:,:)
  real(rt)        , pointer:: q2(:,:,:,:)
  real(rt)        , pointer:: q3(:,:,:,:)


  real(rt)        , pointer:: shk(:,:,:)

  ! temporary interface values of the parabola
  real(rt)        , pointer :: sxm(:,:,:), sym(:,:,:), szm(:,:,:)
  real(rt)        , pointer :: sxp(:,:,:), syp(:,:,:), szp(:,:,:)

  real(rt)        , pointer :: qxm(:,:,:,:), qym(:,:,:,:), qzm(:,:,:,:)
  real(rt)        , pointer :: qxp(:,:,:,:), qyp(:,:,:,:), qzp(:,:,:,:)

  integer :: ngf
  integer :: It_lo(3), It_hi(3)
  integer :: st_lo(3), st_hi(3)
  integer :: shk_lo(3), shk_hi(3)

  real(rt) :: div1
  integer :: i, j, k, n

  type (eos_t) :: eos_state

  ngf = 1

  It_lo = lo(:) - dg(:)
  It_hi = hi(:) + dg(:)

  st_lo = lo(:) - 2*dg(:)
  st_hi = hi(:) + 2*dg(:)

  shk_lo(:) = lo(:) - dg(:)
  shk_hi(:) = hi(:) + dg(:)

  call bl_allocate(   div, lo, hi+dg)

  call bl_allocate(q1, flux1_lo, flux1_hi, NGDNV)
  call bl_allocate(q2, flux2_lo, flux2_hi, NGDNV)


  call bl_allocate(sxm, st_lo, st_hi)
  call bl_allocate(sxp, st_lo, st_hi)
  call bl_allocate(qxm, It_lo, It_hi, NQ)
  call bl_allocate(qxp, It_lo, It_hi, NQ)

  call bl_allocate(sym, st_lo, st_hi)
  call bl_allocate(syp, st_lo, st_hi)
  call bl_allocate(qym, It_lo, It_hi, NQ)
  call bl_allocate(qyp, It_lo, It_hi, NQ)

  call bl_allocate(shk, shk_lo, shk_hi)

  if (ppm_type == 0) then
     call amrex_error("ERROR: method of lines integration does not support ppm_type = 0")
  endif

    ! multidimensional shock detection -- this will be used to do the
    ! hybrid Riemann solver
    if (hybrid_riemann == 1) then
       call shock(q, q_lo, q_hi, shk, shk_lo, shk_hi, lo, hi, dx)
    else
       shk(:,:,:) = ZERO
    endif

  ! Compute flattening coefficient for slope calculations.
  call bl_allocate(flatn, q_lo, q_hi)

  if (first_order_hydro == 1) then
     flatn = ZERO
  elseif (use_flattening == 1) then
     call uflatten(lo - ngf*dg, hi + ngf*dg, &
                   q, flatn, q_lo, q_hi, QPRES)
  else
     flatn = ONE
  endif


  do n = 1, NQ
     call ppm_reconstruct(q, q_lo, q_hi, NQ, n, &
                          flatn, q_lo, q_hi, &
                          sxm, sxp, &
                          sym, syp, &
                          st_lo, st_hi, &
                          lo, hi, dx)

     ! Construct the interface states -- this is essentially just a
     ! reshuffling of interface states from zone-center indexing to
     ! edge-centered indexing
     do k = lo(3)-dg(3), hi(3)+dg(3)
        do j = lo(2)-dg(2), hi(2)+dg(2)
           do i = lo(1)-1, hi(1)+1

              ! x-edges

              ! left state at i-1/2 interface
              qxm(i,j,k,n) = sxp(i-1,j,k)

              ! right state at i-1/2 interface
              qxp(i,j,k,n) = sxm(i,j,k)

              ! y-edges

              ! left state at j-1/2 interface
              qym(i,j,k,n) = syp(i,j-1,k)

              ! right state at j-1/2 interface
              qyp(i,j,k,n) = sym(i,j,k)


           enddo
        enddo
     enddo

     ! use T to define p
     if (ppm_temp_fix == 1) then
        do k = lo(3)-dg(3), hi(3)+dg(3)
           do j = lo(2)-dg(2), hi(2)+dg(2)
              do i = lo(1)-1, hi(1)+1

                 eos_state%rho    = qxp(i,j,k,QRHO)
                 eos_state%T      = qxp(i,j,k,QTEMP)
                 eos_state%xn(:)  = qxp(i,j,k,QFS:QFS-1+nspec)
                 eos_state%aux(:) = qxp(i,j,k,QFX:QFX-1+naux)

                 call eos(eos_input_rt, eos_state)

                 qxp(i,j,k,QPRES) = eos_state%p
                 qxp(i,j,k,QREINT) = qxp(i,j,k,QRHO)*eos_state%e
                 ! should we try to do something about Gamma_! on interface?

                 eos_state%rho    = qxm(i,j,k,QRHO)
                 eos_state%T      = qxm(i,j,k,QTEMP)
                 eos_state%xn(:)  = qxm(i,j,k,QFS:QFS-1+nspec)
                 eos_state%aux(:) = qxm(i,j,k,QFX:QFX-1+naux)

                 call eos(eos_input_rt, eos_state)

                 qxm(i,j,k,QPRES) = eos_state%p
                 qxm(i,j,k,QREINT) = qxm(i,j,k,QRHO)*eos_state%e
                 ! should we try to do something about Gamma_! on interface?

                 eos_state%rho    = qyp(i,j,k,QRHO)
                 eos_state%T      = qyp(i,j,k,QTEMP)
                 eos_state%xn(:)  = qyp(i,j,k,QFS:QFS-1+nspec)
                 eos_state%aux(:) = qyp(i,j,k,QFX:QFX-1+naux)

                 call eos(eos_input_rt, eos_state)

                 qyp(i,j,k,QPRES) = eos_state%p
                 qyp(i,j,k,QREINT) = qyp(i,j,k,QRHO)*eos_state%e
                 ! should we try to do something about Gamma_! on interface?

                 eos_state%rho    = qym(i,j,k,QRHO)
                 eos_state%T      = qym(i,j,k,QTEMP)
                 eos_state%xn(:)  = qym(i,j,k,QFS:QFS-1+nspec)
                 eos_state%aux(:) = qym(i,j,k,QFX:QFX-1+naux)

                 call eos(eos_input_rt, eos_state)

                 qym(i,j,k,QPRES) = eos_state%p
                 qym(i,j,k,QREINT) = qym(i,j,k,QRHO)*eos_state%e
                 ! should we try to do something about Gamma_! on interface?


              end do
           end do
        end do
     endif

     ! Compute F^x at kc (k3d)
     call cmpflx(qxm, qxp, It_lo, It_hi, &
                 flux1, flux1_lo, flux1_hi, &
                 q1, flux1_lo, flux1_hi, &  ! temporary
                 qaux, qa_lo, qa_hi, &
                 shk, shk_lo, shk_hi, &
                 1, [lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)], domlo, domhi)


     ! Compute F^y at kc (k3d)
     call cmpflx(qym, qyp, It_lo, It_hi, &
                 flux2, flux2_lo, flux2_hi, &
                 q2, flux2_lo, flux2_hi, &  ! temporary
                 qaux, qa_lo, qa_hi, &
                 shk, shk_lo, shk_hi, &
                 2, [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)], domlo, domhi)



  end do

  call bl_deallocate(flatn)

  call bl_deallocate(sxm)
  call bl_deallocate(sxp)
  call bl_deallocate(qxm)
  call bl_deallocate(qxp)

  call bl_deallocate(sym)
  call bl_deallocate(syp)
  call bl_deallocate(qym)
  call bl_deallocate(qyp)


  call bl_deallocate(shk)


  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo, hi+dg, q, q_lo, q_hi, &
            dx, div, lo, hi+dg)

  do n = 1, NVAR

     if ( n == UTEMP ) then
        flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
        flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO


     end if

  end do

  call apply_av(flux1_lo, flux1_hi, 1, dx, &
                div, lo, hi+dg, &
                uin, uin_lo, uin_hi, &
                flux1, flux1_lo, flux1_hi)

  call apply_av(flux2_lo, flux2_hi, 2, dx, &
                div, lo, hi+dg, &
                uin, uin_lo, uin_hi, &
                flux2, flux2_lo, flux2_hi)


  if (limit_fluxes_on_small_dens == 1) then
     call limit_hydro_fluxes_on_small_dens(uin,uin_lo,uin_hi, &
                                           q,q_lo,q_hi, &
                                           vol,vol_lo,vol_hi, &
                                           flux1,flux1_lo,flux1_hi, &
                                           area1,area1_lo,area1_hi, &
                                           flux2,flux2_lo,flux2_hi, &
                                           area2,area2_lo,area2_hi, &
                                           lo,hi,dt,dx)

  endif

  call normalize_species_fluxes(flux1_lo, flux1_hi, flux1, flux1_lo, flux1_hi)
  call normalize_species_fluxes(flux2_lo, flux2_hi, flux2, flux2_lo, flux2_hi)


  ! For hydro, we will create an update source term that is
  ! essentially the flux divergence.  This can be added with dt to
  ! get the update
  do n = 1, NVAR
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)

              update(i,j,k,n) = update(i,j,k,n) + &
                   (flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) + &
                    flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) ) / vol(i,j,k)



              if (n == UMX) then
                 ! add the pressure source term for axisymmetry
                 if (coord_type > 0) then
                    update(i,j,k,n) = update(i,j,k,n) - (q1(i+1,j,k,GDPRES) - q1(i,j,k,GDPRES))/ dx(1)
                 endif
              endif

              ! for storage
              update_flux(i,j,k,n) = update_flux(i,j,k,n) + &
                   stage_weight * update(i,j,k,n)

              update(i,j,k,n) = update(i,j,k,n) + srcU(i,j,k,n)

           enddo
        enddo
     enddo
  enddo




  ! Scale the fluxes for the form we expect later in refluxing.

  call scale_flux(flux1_lo, flux1_hi, flux1, flux1_lo, flux1_hi, area1, area1_lo, area1_hi, dt)
  call scale_flux(flux2_lo, flux2_hi, flux2, flux2_lo, flux2_hi, area2, area2_lo, area2_hi, dt)


  if (coord_type > 0) then
     pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = q1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),GDPRES) * dt
  end if

  call bl_deallocate(   div)

  call bl_deallocate(q1)
  call bl_deallocate(q2)


end subroutine ca_mol_single_stage



module mol_module_cuda

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  subroutine ca_construct_flux_cuda(lo, hi, domlo, domhi, dx, dt, idir, &
                                    uin, uin_lo, uin_hi, &
                                    div, div_lo, div_hi, &
                                    qaux, qa_lo, qa_hi, &
                                    qm, qm_lo, qm_hi, &
                                    qp, qp_lo, qp_hi, &
                                    qint, qe_lo, qe_hi, &
                                    flux, f_lo, f_hi, &
                                    area, a_lo, a_hi) &
                                    bind(c,name='ca_construct_flux_cuda')

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, NGDNV, NQAUX, NQ
    use advection_util_module, only: apply_av, normalize_species_fluxes, scale_flux
    use riemann_module, only: cmpflx_cuda

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ), value :: idir
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: uin_lo(3), uin_hi(3)
    integer,  intent(in   ) :: div_lo(3), div_hi(3)
    integer,  intent(in   ) :: qa_lo(3), qa_hi(3)
    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)
    integer,  intent(in   ) :: qe_lo(3), qe_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: a_lo(3), a_hi(3)

    real(rt), intent(in   ) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt), intent(in   ) :: div(div_lo(1):div_hi(1), div_lo(2):div_hi(2), div_lo(3):div_hi(3))
    real(rt), intent(in   ) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ,3)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ,3)
    real(rt), intent(inout) :: qint(qe_lo(1):qe_hi(1), qe_lo(2):qe_hi(2), qe_lo(3):qe_hi(3), NGDNV)
    real(rt), intent(inout) :: flux(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3), NVAR)
    real(rt), intent(in   ) :: area(a_lo(1):a_hi(1), a_lo(2):a_hi(2), a_lo(3):a_hi(3))
    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in   ), value :: dt

    !$gpu

    call cmpflx_cuda(lo, hi, domlo, domhi, idir, qm, qm_lo, qm_hi, qp, qp_lo, qp_hi, &
                     qint, qe_lo, qe_hi, flux, f_lo, f_hi, qaux, qa_lo, qa_hi)
    call apply_av(lo, hi, idir, dx, div, div_lo, div_hi, uin, uin_lo, uin_hi, flux, f_lo, f_hi)
    call normalize_species_fluxes(lo, hi, flux, f_lo, f_hi)
    call scale_flux(lo, hi, flux, f_lo, f_hi, area, a_lo, a_hi, dt)

  end subroutine ca_construct_flux_cuda

end module mol_module_cuda
