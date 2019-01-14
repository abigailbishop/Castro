












! advection routines in support of the CTU unsplit advection scheme

module ctu_module

  implicit none

contains

  subroutine consup(uin, uin_lo, uin_hi, &
                    q, q_lo, q_hi, &
                    uout, uout_lo, uout_hi, &
                    update, updt_lo, updt_hi, &
                    flux1, flux1_lo, flux1_hi, &
                    flux2, flux2_lo, flux2_hi, &
                    qx, qx_lo, qx_hi, &
                    qy, qy_lo, qy_hi, &
                    area1, area1_lo, area1_hi, &
                    area2, area2_lo, area2_hi, &
                    vol,vol_lo,vol_hi, &
                    div, lo, hi, dx, dt, &
                    mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                    eden_lost, xang_lost, yang_lost, zang_lost, &
                    verbose)

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, NGDNV, NQ, &
                                   GDPRES, &
                                   track_grid_losses, limit_fluxes_on_small_dens
    use advection_util_module, only : limit_hydro_fluxes_on_small_dens, normalize_species_fluxes, calc_pdivu
    use castro_util_module, only : position, linear_to_angular_momentum
    use prob_params_module, only : mom_flux_has_p, domlo_level, domhi_level, center, dg, coord_type
    use amrinfo_module, only : amr_level
    use amrex_constants_module, only : ZERO, ONE, TWO, FOURTH, HALF

    use amrex_fort_module, only : rt => amrex_real

    integer, intent(in) ::       lo(3),       hi(3)
    integer, intent(in) ::   uin_lo(3),   uin_hi(3)
    integer, intent(in) ::     q_lo(3),     q_hi(3)
    integer, intent(in) ::  uout_lo(3),  uout_hi(3)
    integer, intent(in) ::  updt_lo(3),  updt_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) ::    qy_lo(3),    qy_hi(3)
    integer, intent(in) ::    qx_lo(3),    qx_hi(3)
    integer, intent(in) ::   vol_lo(3),   vol_hi(3)


    integer, intent(in) :: verbose


    real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1),uin_lo(2):uin_hi(2),uin_lo(3):uin_hi(3),NVAR)
    real(rt)        , intent(in) :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
    real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1),updt_lo(2):updt_hi(2),updt_lo(3):updt_hi(3),NVAR)

    real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1),flux1_lo(2):flux1_hi(2),flux1_lo(3):flux1_hi(3),NVAR)
    real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1),area1_lo(2):area1_hi(2),area1_lo(3):area1_hi(3))
    real(rt)        , intent(in) ::    qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)

    real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1),flux2_lo(2):flux2_hi(2),flux2_lo(3):flux2_hi(3),NVAR)
    real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1),area2_lo(2):area2_hi(2),area2_lo(3):area2_hi(3))
    real(rt)        , intent(in) ::    qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)


    real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt)        , intent(in) :: div(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
    real(rt)        , intent(in) :: dx(3), dt


    real(rt)        , intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
    real(rt)        , intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

    real(rt)         :: div1, volinv
    integer          :: i, j, g, k, n
    integer          :: domlo(3), domhi(3)
    real(rt)         :: loc(3), ang_mom(3)

    real(rt)        , pointer:: pdivu(:,:,:)


    call bl_allocate(pdivu, lo, hi)

    call calc_pdivu(lo, hi, &
                    qx, qx_lo, qx_hi, &
                    area1, area1_lo, area1_hi, &
                    qy, qy_lo, qy_hi, &
                    area2, area2_lo, area2_hi, &
                    vol, vol_lo, vol_hi, &
                    dx, pdivu, lo, hi)

    do n = 1, NVAR

       if ( n == UTEMP ) then

          flux1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),n) = ZERO
          flux2(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3),n) = ZERO


       else

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)+1
                   div1 = FOURTH*(div(i,j,k) + div(i,j+dg(2),k) + &
                                  div(i,j,k+dg(3)) + div(i,j+dg(2),k+dg(3)))
                   div1 = difmag*min(ZERO, div1)

                   flux1(i,j,k,n) = flux1(i,j,k,n) + dx(1) * div1 * (uin(i,j,k,n)-uin(i-1,j,k,n))
                enddo
             enddo
          enddo

          do k = lo(3),hi(3)
             do j = lo(2),hi(2)+1
                do i = lo(1),hi(1)
                   div1 = FOURTH*(div(i,j,k) + div(i+1,j,k) + &
                                  div(i,j,k+dg(3)) + div(i+1,j,k+dg(3)))
                   div1 = difmag*min(ZERO, div1)

                   flux2(i,j,k,n) = flux2(i,j,k,n) + dx(2) * div1 * (uin(i,j,k,n)-uin(i,j-1,k,n))
                enddo
             enddo
          enddo


       endif

    enddo


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

    call normalize_species_fluxes(flux1_lo, flux1_hi, flux1, flux1_lo,flux1_hi)
    call normalize_species_fluxes(flux2_lo, flux2_hi, flux2, flux2_lo,flux2_hi)

    ! For hydro, we will create an update source term that is
    ! essentially the flux divergence.  This can be added with dt to
    ! get the update

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                volinv = ONE / vol(i,j,k)

                update(i,j,k,n) = update(i,j,k,n) + &
                     ( flux1(i,j,k,n) * area1(i,j,k) - flux1(i+1,j,k,n) * area1(i+1,j,k) &
                     + flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j+1,k,n) * area2(i,j+1,k) &
                     ) * volinv

                ! Add the p div(u) source term to (rho e).
                if (n .eq. UEINT) then
                   update(i,j,k,n) = update(i,j,k,n) - pdivu(i,j,k)
                endif

             enddo
          enddo
       enddo
    enddo



    ! Add gradp term to momentum equation -- only for axisymmetric
    ! coords (and only for the radial flux).

    if (.not. mom_flux_has_p(1)%comp(UMX)) then
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                update(i,j,k,UMX) = update(i,j,k,UMX) - (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx(1)
                !update(i,j,UMY) = update(i,j,UMY) - (pgdy(i,j+1)-pgdy(i,j)) / dy
             enddo
          enddo
       enddo
    endif


    ! Scale the fluxes for the form we expect later in refluxing.

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1) + 1
                flux1(i,j,k,n) = dt * flux1(i,j,k,n) * area1(i,j,k)
             enddo
          enddo
       enddo
    enddo

    do n = 1, NVAR
       do k = lo(3), hi(3)
          do j = lo(2), hi(2) + 1
             do i = lo(1), hi(1)
                flux2(i,j,k,n) = dt * flux2(i,j,k,n) * area2(i,j,k)
             enddo
          enddo
       enddo
    enddo




    ! Add up some diagnostic quantities. Note that we are not dividing by the cell volume.

    if (track_grid_losses .eq. 1) then

       domlo = domlo_level(:,amr_level)
       domhi = domhi_level(:,amr_level)


       if (lo(2) .le. domlo(2) .and. hi(2) .ge. domlo(2)) then

          j = domlo(2)
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccy=.false.)

                mass_lost = mass_lost - flux2(i,j,k,URHO)
                xmom_lost = xmom_lost - flux2(i,j,k,UMX)
                ymom_lost = ymom_lost - flux2(i,j,k,UMY)
                zmom_lost = zmom_lost - flux2(i,j,k,UMZ)
                eden_lost = eden_lost - flux2(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux2(i,j,k,UMX:UMZ))
                xang_lost = xang_lost - ang_mom(1)
                yang_lost = yang_lost - ang_mom(2)
                zang_lost = zang_lost - ang_mom(3)

             enddo
          enddo

       endif

       if (lo(2) .le. domhi(2) .and. hi(2) .ge. domhi(2)) then

          j = domhi(2) + 1
          do k = lo(3), hi(3)
             do i = lo(1), hi(1)

                loc = position(i,j,k,ccy=.false.)

                mass_lost = mass_lost + flux2(i,j,k,URHO)
                xmom_lost = xmom_lost + flux2(i,j,k,UMX)
                ymom_lost = ymom_lost + flux2(i,j,k,UMY)
                zmom_lost = zmom_lost + flux2(i,j,k,UMZ)
                eden_lost = eden_lost + flux2(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux2(i,j,k,UMX:UMZ))
                xang_lost = xang_lost + ang_mom(1)
                yang_lost = yang_lost + ang_mom(2)
                zang_lost = zang_lost + ang_mom(3)

             enddo
          enddo

       endif

       if (lo(1) .le. domlo(1) .and. hi(1) .ge. domlo(1)) then

          i = domlo(1)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                loc = position(i,j,k,ccx=.false.)

                mass_lost = mass_lost - flux1(i,j,k,URHO)
                xmom_lost = xmom_lost - flux1(i,j,k,UMX)
                ymom_lost = ymom_lost - flux1(i,j,k,UMY)
                zmom_lost = zmom_lost - flux1(i,j,k,UMZ)
                eden_lost = eden_lost - flux1(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux1(i,j,k,UMX:UMZ))
                xang_lost = xang_lost - ang_mom(1)
                yang_lost = yang_lost - ang_mom(2)
                zang_lost = zang_lost - ang_mom(3)

             enddo
          enddo

       endif

       if (lo(1) .le. domhi(1) .and. hi(1) .ge. domhi(1)) then

          i = domhi(1) + 1
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)

                loc = position(i,j,k,ccx=.false.)

                mass_lost = mass_lost + flux1(i,j,k,URHO)
                xmom_lost = xmom_lost + flux1(i,j,k,UMX)
                ymom_lost = ymom_lost + flux1(i,j,k,UMY)
                zmom_lost = zmom_lost + flux1(i,j,k,UMZ)
                eden_lost = eden_lost + flux1(i,j,k,UEDEN)

                ang_mom   = linear_to_angular_momentum(loc - center, flux1(i,j,k,UMX:UMZ))
                xang_lost = xang_lost + ang_mom(1)
                yang_lost = yang_lost + ang_mom(2)
                zang_lost = zang_lost + ang_mom(3)

             enddo
          enddo

       endif

    endif

    call bl_deallocate(pdivu)

  end subroutine consup


  subroutine ca_ctu_update(lo, hi, is_finest_level, time, &
                           domlo, domhi, &
                           uin, uin_lo, uin_hi, &
                           uout, uout_lo, uout_hi, &
                           q, q_lo, q_hi, &
                           qaux, qa_lo, qa_hi, &
                           srcQ, srQ_lo, srQ_hi, &
                           update, updt_lo, updt_hi, &
                           delta, dt, &
                           flux1, flux1_lo, flux1_hi, &
                           flux2, flux2_lo, flux2_hi, &
                           area1, area1_lo, area1_hi, &
                           area2, area2_lo, area2_hi, &
                           pradial, p_lo, p_hi, &
                           dloga, dloga_lo, dloga_hi, &
                           vol, vol_lo, vol_hi, &
                           verbose, &
                           mass_lost, xmom_lost, ymom_lost, zmom_lost, &
                           eden_lost, xang_lost, yang_lost, zang_lost) bind(C, name="ca_ctu_update")

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : NQ, QVAR, QPRES, NQAUX, NVAR, NHYP, NGDNV, UMX, GDPRES, &
                                   use_flattening, &
                                   first_order_hydro
    use advection_util_module, only : divu
    use amrex_constants_module, only : ZERO, ONE
    use flatten_module, only: uflatten
    use prob_params_module, only : mom_flux_has_p, dg, coord_type
    use ctu_advection_module, only : umeth

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: is_finest_level
    integer, intent(in) :: lo(3), hi(3), verbose
    integer, intent(in) :: domlo(3), domhi(3)
    integer, intent(in) :: uin_lo(3), uin_hi(3)
    integer, intent(in) :: uout_lo(3), uout_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: srQ_lo(3), srQ_hi(3)
    integer, intent(in) :: updt_lo(3), updt_hi(3)
    integer, intent(in) :: flux1_lo(3), flux1_hi(3)
    integer, intent(in) :: flux2_lo(3), flux2_hi(3)
    integer, intent(in) :: area1_lo(3), area1_hi(3)
    integer, intent(in) :: area2_lo(3), area2_hi(3)
    integer, intent(in) :: vol_lo(3), vol_hi(3)
    integer, intent(in) :: p_lo(3), p_hi(3)
    integer, intent(in) :: dloga_lo(3), dloga_hi(3)

    real(rt)        , intent(in) :: uin(uin_lo(1):uin_hi(1), uin_lo(2):uin_hi(2), uin_lo(3):uin_hi(3), NVAR)
    real(rt)        , intent(inout) :: uout(uout_lo(1):uout_hi(1), uout_lo(2):uout_hi(2), uout_lo(3):uout_hi(3), NVAR)
    real(rt)        , intent(inout) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt)        , intent(in) :: qaux(qa_lo(1):qa_hi(1), qa_lo(2):qa_hi(2), qa_lo(3):qa_hi(3), NQAUX)
    real(rt)        , intent(in) :: srcQ(srQ_lo(1):srQ_hi(1), srQ_lo(2):srQ_hi(2), srQ_lo(3):srQ_hi(3), QVAR)
    real(rt)        , intent(inout) :: update(updt_lo(1):updt_hi(1), updt_lo(2):updt_hi(2), updt_lo(3):updt_hi(3), NVAR)
    real(rt)        , intent(inout) :: flux1(flux1_lo(1):flux1_hi(1), flux1_lo(2):flux1_hi(2), flux1_lo(3):flux1_hi(3), NVAR)
    real(rt)        , intent(inout) :: flux2(flux2_lo(1):flux2_hi(1), flux2_lo(2):flux2_hi(2), flux2_lo(3):flux2_hi(3), NVAR)
    real(rt)        , intent(in) :: area1(area1_lo(1):area1_hi(1), area1_lo(2):area1_hi(2), area1_lo(3):area1_hi(3))
    real(rt)        , intent(in) :: area2(area2_lo(1):area2_hi(1), area2_lo(2):area2_hi(2), area2_lo(3):area2_hi(3))
    real(rt)        , intent(in) :: vol(vol_lo(1):vol_hi(1), vol_lo(2):vol_hi(2), vol_lo(3):vol_hi(3))

    real(rt)        , intent(in) :: dloga(dloga_lo(1):dloga_hi(1),dloga_lo(2):dloga_hi(2),dloga_lo(3):dloga_hi(3))
    real(rt)        , intent(inout) :: pradial(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3))

    real(rt)        , intent(in) :: delta(3), dt, time

    real(rt)        , intent(inout) :: mass_lost, xmom_lost, ymom_lost, zmom_lost
    real(rt)        , intent(inout) :: eden_lost, xang_lost, yang_lost, zang_lost

    ! Automatic arrays for workspace
    real(rt)        , pointer:: flatn(:,:,:)
    real(rt)        , pointer:: div(:,:,:)

    ! Edge-centered primitive variables (Riemann state)
    real(rt)        , pointer:: q1(:,:,:,:)
    real(rt)        , pointer:: q2(:,:,:,:)
    real(rt)        , pointer:: q3(:,:,:,:)

    integer :: ngf
    integer :: q1_lo(3), q1_hi(3), q2_lo(3), q2_hi(3), q3_lo(3), q3_hi(3)

    ngf = 1

    call bl_allocate(   div, lo, hi+dg)

    q1_lo = flux1_lo - dg
    q1_hi = flux1_hi + dg
    q2_lo = flux2_lo - dg
    q2_hi = flux2_hi + dg

    call bl_allocate(q1, q1_lo, q1_hi, NGDNV)
    call bl_allocate(q2, q2_lo, q2_hi, NGDNV)

    ! Compute flattening coefficient for slope calculations.
    call bl_allocate( flatn, q_lo, q_hi)

    if (first_order_hydro == 1) then
       flatn = ZERO
    elseif (use_flattening == 1) then
       call uflatten(lo-dg*ngf, hi+dg*ngf, &
                     q, flatn, q_lo, q_hi, QPRES)
    else
       flatn = ONE
    endif

    ! Compute hyperbolic fluxes using unsplit Godunov
    call umeth(q, q_lo, q_hi, &
               flatn, &
               qaux, qa_lo, qa_hi, &
               srcQ, srQ_lo, srQ_hi, &
               lo, hi, delta, dt, &
               uout, uout_lo, uout_hi, &
               flux1, flux1_lo, flux1_hi, &
               flux2, flux2_lo, flux2_hi, &
               q1, q1_lo, q1_hi, &
               q2, q2_lo, q2_hi, &
                area1, area1_lo, area1_hi, &
                area2, area2_lo, area2_hi, &
                vol, vol_lo, vol_hi, &
                dloga, dloga_lo, dloga_hi, &
                domlo, domhi)


    call bl_deallocate( flatn)

    ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
    call divu(lo, hi+dg, q, q_lo, q_hi, delta, div, lo, hi+dg)

    ! Conservative update
    call consup(uin, uin_lo, uin_hi, &
                q, q_lo, q_hi, &
                uout, uout_lo, uout_hi, &
                update, updt_lo, updt_hi, &
                flux1, flux1_lo, flux1_hi, &
                flux2, flux2_lo, flux2_hi, &
                q1, q1_lo, q1_hi, &
                q2, q2_lo, q2_hi, &
                area1, area1_lo, area1_hi, &
                area2, area2_lo, area2_hi, &
                vol, vol_lo, vol_hi, &
                div, lo, hi, delta, dt, &
                mass_lost,xmom_lost,ymom_lost,zmom_lost, &
                eden_lost,xang_lost,yang_lost,zang_lost, &
                verbose)


    if (.not. mom_flux_has_p(1)%comp(UMX)) then
       pradial(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = q1(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3),GDPRES) * dt
    end if

    call bl_deallocate(   div)

    call bl_deallocate(    q1)
    call bl_deallocate(    q2)

  end subroutine ca_ctu_update

end module ctu_module
