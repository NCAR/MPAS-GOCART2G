! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_microphysics
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer, only : mpas_timer_start, mpas_timer_stop

 use mpas_atmphys_constants
 use mpas_atmphys_init_microphysics
 use mpas_atmphys_interface
 use mpas_atmphys_vars

#ifdef DO_GOCART2G
 use mpas_chemistry_gocart2G_vars
#endif

!wrf physics:
 use module_mp_kessler
 use module_mp_thompson
 use module_mp_wsm6,only: wsm6
 use mp_wsm6,only: mp_wsm6_init,refl10cm_wsm6


 implicit none
 private
 public:: allocate_microphysics,   &
          deallocate_microphysics, &
          driver_microphysics,     &
          init_microphysics


!MPAS driver for parameterization of cloud microphysics processes.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_driver_microphysics:
! ------------------------------------------------
! allocate_microphysics     : allocate local arrays for parameterization of cloud microphysics.
! deallocate_microphysics   : deallocate local arrays for parameterization of cloud microphysics.
! microphysics_init         : initialization of individual cloud microphysics schemes.
! driver_microphysics       : main driver (called from mpas_atm_time_integration).
! precip_from_MPAS          : initialize timestep local arrays for precipitation.
! precip_to_MPAS            : copy local arrays to MPAS arrays.
! compute_radar_reflectivity: compute radar reflectivities.
! compute_relhum            : compute relative humidity.
!
! WRF physics called from microphysics_driver:
! --------------------------------------------
! * module_mp_kessler : Kessler cloud microphysics.
! * module_mp_thompson: Thompson cloud microphysics.
! * module_mp_wsm6    : WSM6 cloud microphysics.
!
! comments:
! unlike all the other physics parameterizations, parameterizations of cloud microphysics schemes
! are called at the bottom of subroutine atm_srk3 in module atm_time_integration.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed call to the Thompson cloud microphysics scheme until scheme is updated to that in WRF revision 3.5.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * added subroutine compute_relhum to calculate the relative humidity using the functions rslf and rsif from
!   the Thompson cloud microphysics scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-12.
! * removed the argument tend from the call to microphysics_from_MPAS (not needed).
!   Laura D. Fowler (laura@ucar.edu) / 2013-11-07.
! * in call to subroutine wsm6, replaced the variable g (that originally pointed to gravity) with gravity,
!   for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * moved the variable relhum from the diag_physics to the diag pool. Changed the argument
!   list for the subroutine compute_relhum accordingly.
!   Laura D. Fowler (laura@ucar.edu) / 2015-04-22.
! * added parameterization of the Thompson cloud microphysics from WRF version 3.8.
!   Laura D. Fowler (laura@ucar.edu) / 2016-03-28.
! * in subroutine compute_relhum, multiply relhum by 100. so that it has the same unit as in the initial
!   conditions.
!   Laura D. Fowler (laura@ucar.edu) / 2016-06-20.
! * added parameterization of the WSM6 cloud microphysics from WRF version 3.8.1. To initialize WSM6 as in its
!   original version, set the hail_option to 0.
!   Laura D. Fowler (laura@ucar.edu) / 2016-09-19.
! * since we removed the local variable microp_scheme from mpas_atmphys_vars.F, now defines microp_scheme as a
!   pointer to config_microp_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2917-02-16.

!--- initialization option for WSM6 from WRF version 3.8.1. this option could also be set as a namelist parameter.
 integer,parameter:: hail_opt = 0


 contains


!=================================================================================================================
 subroutine allocate_microphysics(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

!sounding variables:
 if(.not.allocated(rho_p) ) allocate(rho_p(ims:ime,kms:kme,jms:jme) )
 if(.not.allocated(th_p)  ) allocate(th_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pi_p)  ) allocate(pi_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(pres_p)) allocate(pres_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(z_p)   ) allocate(z_p(ims:ime,kms:kme,jms:jme)   )
 if(.not.allocated(dz_p)  ) allocate(dz_p(ims:ime,kms:kme,jms:jme)  )
 if(.not.allocated(w_p)   ) allocate(w_p(ims:ime,kms:kme,jms:jme)   )

!mass mixing ratios:
 if(.not.allocated(qv_p)) allocate(qv_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(qc_p)) allocate(qc_p(ims:ime,kms:kme,jms:jme))
 if(.not.allocated(qr_p)) allocate(qr_p(ims:ime,kms:kme,jms:jme))

 !surface precipitation:
 if(.not.allocated(rainnc_p) ) allocate(rainnc_p(ims:ime,jms:jme) )
 if(.not.allocated(rainncv_p)) allocate(rainncv_p(ims:ime,jms:jme))

 microp_select: select case(trim(microp_scheme))
    case ("mp_thompson","mp_thompson_aerosols","mp_thompson_gocart2G","mp_wsm6")
       !mass mixing ratios:
       if(.not.allocated(qi_p)) allocate(qi_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(qs_p)) allocate(qs_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(qg_p)) allocate(qg_p(ims:ime,kms:kme,jms:jme))

       !surface precipitation:
       if(.not.allocated(sr_p)        ) allocate(sr_p(ims:ime,jms:jme)        )
       if(.not.allocated(snownc_p)    ) allocate(snownc_p(ims:ime,jms:jme)    )
       if(.not.allocated(snowncv_p)   ) allocate(snowncv_p(ims:ime,jms:jme)   )
       if(.not.allocated(graupelnc_p) ) allocate(graupelnc_p(ims:ime,jms:jme) )
       if(.not.allocated(graupelncv_p)) allocate(graupelncv_p(ims:ime,jms:jme))

       !cloud water,cloud ice,and snow effective radii:
       if(.not.allocated(recloud_p)) allocate(recloud_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(reice_p)  ) allocate(reice_p(ims:ime,kms:kme,jms:jme)  )
       if(.not.allocated(resnow_p) ) allocate(resnow_p(ims:ime,kms:kme,jms:jme) )

       !precipitation flux:
       if(.not.allocated(rainprod_p)) allocate(rainprod_p(ims:ime,kms:kme,jms:jme))
       if(.not.allocated(evapprod_p)) allocate(evapprod_p(ims:ime,kms:kme,jms:jme))

    microp2_select: select case(trim(microp_scheme))
       case("mp_thompson","mp_thompson_aerosols","mp_thompson_gocart2G")
          if(.not.allocated(ntc_p)) allocate(ntc_p(ims:ime,jms:jme))
          if(.not.allocated(muc_p)) allocate(muc_p(ims:ime,jms:jme))
          if(.not.allocated(ni_p) ) allocate(ni_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(nr_p) ) allocate(nr_p(ims:ime,kms:kme,jms:jme))

          if(.not.allocated(prcpfluxr_p)) allocate(prcpfluxr_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(prcpfluxi_p)) allocate(prcpfluxi_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(prcpfluxs_p)) allocate(prcpfluxs_p(ims:ime,kms:kme,jms:jme))
          if(.not.allocated(prcpfluxg_p)) allocate(prcpfluxg_p(ims:ime,kms:kme,jms:jme))

          microp3_select: select case(trim(microp_scheme))
             case("mp_thompson_aerosols")
                if(.not.allocated(nifa2d_p)) allocate(nifa2d_p(ims:ime,jms:jme))
                if(.not.allocated(nwfa2d_p)) allocate(nwfa2d_p(ims:ime,jms:jme))
                if(.not.allocated(nc_p)    ) allocate(nc_p(ims:ime,kms:kme,jms:jme)  )
                if(.not.allocated(nifa_p)  ) allocate(nifa_p(ims:ime,kms:kme,jms:jme))
                if(.not.allocated(nwfa_p)  ) allocate(nwfa_p(ims:ime,kms:kme,jms:jme))

              case("mp_thompson_gocart2G")
                if(.not.allocated(nc_p)  ) allocate(nc_p(ims:ime,kms:kme,jms:jme)  )
                if(.not.allocated(nifa_p)) allocate(nifa_p(ims:ime,kms:kme,jms:jme))
                if(.not.allocated(nwfa_p)) allocate(nwfa_p(ims:ime,kms:kme,jms:jme))

             case default
          end select microp3_select

       case default
    end select microp2_select

    case default
 end select microp_select

 end subroutine allocate_microphysics

!=================================================================================================================
 subroutine deallocate_microphysics(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

!sounding variables:
 if(allocated(rho_p) ) deallocate(rho_p )
 if(allocated(th_p)  ) deallocate(th_p  )
 if(allocated(pi_p)  ) deallocate(pi_p  )
 if(allocated(pres_p)) deallocate(pres_p)
 if(allocated(z_p)   ) deallocate(z_p   )
 if(allocated(dz_p)  ) deallocate(dz_p  )
 if(allocated(w_p)   ) deallocate(w_p   )

!mass mixing ratios:
 if(allocated(qv_p)) deallocate(qv_p)
 if(allocated(qc_p)) deallocate(qc_p)
 if(allocated(qr_p)) deallocate(qr_p)

 !surface precipitation:
 if(allocated(rainnc_p) ) deallocate(rainnc_p )
 if(allocated(rainncv_p)) deallocate(rainncv_p)

 microp_select: select case(trim(microp_scheme))
    case ("mp_thompson","mp_thompson_aerosols","mp_thompson_gocart2G","mp_wsm6")
       !mass mixing ratios:
       if(allocated(qi_p)) deallocate(qi_p)
       if(allocated(qs_p)) deallocate(qs_p)
       if(allocated(qg_p)) deallocate(qg_p)

       !surface precipitation:
       if(allocated(sr_p)        ) deallocate(sr_p        )
       if(allocated(snownc_p)    ) deallocate(snownc_p    )
       if(allocated(snowncv_p)   ) deallocate(snowncv_p   )
       if(allocated(graupelnc_p) ) deallocate(graupelnc_p )
       if(allocated(graupelncv_p)) deallocate(graupelncv_p)

       !cloud water,cloud ice,and snow effective radii:
       if(allocated(recloud_p)) deallocate(recloud_p)
       if(allocated(reice_p)  ) deallocate(reice_p  )
       if(allocated(resnow_p) ) deallocate(resnow_p )

       !precipitation flux:
       if(allocated(rainprod_p)) deallocate(rainprod_p)
       if(allocated(evapprod_p)) deallocate(evapprod_p)

    microp2_select: select case(trim(microp_scheme))
       case("mp_thompson","mp_thompson_aerosols","mp_thompson_gocart2G")
          if(allocated(ntc_p)) deallocate(ntc_p)
          if(allocated(muc_p)) deallocate(muc_p)
          if(allocated(ni_p) ) deallocate(ni_p )
          if(allocated(nr_p) ) deallocate(nr_p )

          if(allocated(prcpfluxr_p)) deallocate(prcpfluxr_p)
          if(allocated(prcpfluxi_p)) deallocate(prcpfluxi_p)
          if(allocated(prcpfluxs_p)) deallocate(prcpfluxs_p)
          if(allocated(prcpfluxg_p)) deallocate(prcpfluxg_p)

          microp3_select: select case(trim(microp_scheme))
             case("mp_thompson_aerosols")
                if(allocated(nifa2d_p)) deallocate(nifa2d_p)
                if(allocated(nwfa2d_p)) deallocate(nwfa2d_p)
                if(allocated(nc_p)    ) deallocate(nc_p    )
                if(allocated(nifa_p)  ) deallocate(nifa_p  )
                if(allocated(nwfa_p)  ) deallocate(nwfa_p  )

              case("mp_thompson_gocart2G")
                if(allocated(nc_p)  ) deallocate(nc_p  )
                if(allocated(nifa_p)) deallocate(nifa_p)
                if(allocated(nwfa_p)) deallocate(nwfa_p)

             case default
          end select microp3_select

       case default
    end select microp2_select

    case default
 end select microp_select

 end subroutine deallocate_microphysics

!=================================================================================================================
 subroutine init_microphysics(dminfo,configs,mesh,state,time_lev,sfc_input,diag_physics)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input
 integer,intent(in):: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: state

!local pointer:
 logical:: do_gocart2G
 logical,pointer:: do_restart
 character(len=StrKIND),pointer:: microp_scheme

!CCPP-compliant flags:
 character(len=StrKIND):: errmsg
 integer:: errflg

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
!call mpas_log_write('--- enter subroutine init_microphysics:')

!initialization of CCPP-compliant flags:
 errmsg = ' '
 errflg = 0

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)
 call mpas_pool_get_config(configs,'config_do_restart'   ,do_restart   )

 microp_select: select case(trim(microp_scheme))
    case("mp_thompson","mp_thompson_aerosols","mp_thompson_gocart2G")
       call thompson_init(l_mp_tables)
       call init_thompson_clouddroplets_forMPAS(mesh,sfc_input,diag_physics)

       microp2_select: select case(trim(microp_scheme))
          case("mp_thompson_aerosols","mp_thompson_gocart2G")
             do_gocart2G = .false.
             if(trim(microp_scheme) == "mp_thompson_gocart2G") do_gocart2G = .true.
             call init_thompson_aerosols_forMPAS(do_restart,do_gocart2G,dminfo,mesh,state,time_lev,diag_physics)

          case default
       end select microp2_select

    case("mp_wsm6")
       call mp_wsm6_init(den0=rho_a,denr=rho_r,dens=rho_s,cl=cliq,cpv=cpv, &
                         hail_opt=hail_opt,errmsg=errmsg,errflg=errflg)

    case default
 end select microp_select

!call mpas_log_write('--- end subroutine init_microphysics:')

 end subroutine init_microphysics

!=================================================================================================================
 subroutine driver_microphysics(configs,mesh,state,time_lev,diag,diag_physics,tend_physics, &
                                tend,itimestep,block_structs,its,ite)
!=================================================================================================================

 use mpas_constants, only : rvord

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: time_lev
 integer,intent(in):: itimestep
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics
 type(mpas_pool_type),intent(inout):: tend

 type(mpas_pool_type), intent(inout) :: block_structs
 type (mpas_pool_type), pointer :: NI2G_diags, SU2G_diags, CA2G_bc_diags, &
                                   CA2G_br_diags, CA2G_oc_diags, DU2G_diags, SS2G_diags
 real(kind=RKIND), dimension(:), pointer :: niWT_bin1, niWT_bin2, niWT_bin3, nh3WT, nh4WT, &
                                            suWT_so2, suWT_so4, suWT_dms, suWT_msa, &
                                            bcWT_phobic, bcWT_philic, brWT_phobic, brWT_philic, &
                                            ocWT_phobic, ocWT_philic, duWT_bin1, duWT_bin2, &
                                            duWT_bin3, duWT_bin4, duWT_bin5, ssWT_bin1, &
                                            ssWT_bin2, ssWT_bin3, ssWT_bin4, ssWT_bin5
 integer :: ic, idx, offset
 real(kind=RKIND) :: col_mass_before, col_mass_after

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme

!local variables and arrays:
 integer:: istep
 integer:: i,j,k

!CCPP-compliant flags:
 character(len=StrKIND):: errmsg
 integer:: errflg

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('---enter subroutine driver_microphysics:')

!initialization of CCPP-compliant flags:
 errmsg = ' '
 errflg = 0

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

!... allocation of microphysics arrays:
!$OMP MASTER
 call allocate_microphysics(configs)
!$OMP END MASTER
!$OMP BARRIER

!... initialization of precipitation related arrays:
 call precip_from_MPAS(configs,diag_physics,its,ite)

!... initialization of soundings for non-hydrostatic dynamical cores.
 call microphysics_from_MPAS(configs,mesh,state,time_lev,diag,diag_physics,tend_physics,its,ite)

!... call to different cloud microphysics schemes:
 microp_select: select case(trim(microp_scheme))
    case ("mp_kessler")
       call mpas_timer_start('mp_kessler')
       call kessler( &
                  t        = th_p      , qv    = qv_p  , qc     = qc_p     ,                &
                  qr       = qr_p      , rho   = rho_p , pii    = pi_p     ,                &
                  dt_in    = dt_microp , z     = z_p   , xlv    = xlv      ,                &
                  cp       = cp        , ep2   = ep_2  , svp1   = svp1     ,                &
                  svp2     = svp2      , svp3  = svp3  , svpt0  = svpt0    ,                &
                  rhowater = rho_w     , dz8w  = dz_p  , rainnc = rainnc_p ,                &
                  rainncv  = rainncv_p ,                                                    &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde   , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme   , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte     &
                   )
       call mpas_timer_stop('mp_kessler')

    case ("mp_thompson")
       call mpas_timer_start('mp_thompson')
       istep = 1
       do while (istep .le. n_microp)
       call mp_gt_driver( &
                  th       = th_p        , qv        = qv_p        , qc         = qc_p         , &
                  qr       = qr_p        , qi        = qi_p        , qs         = qs_p         , &
                  qg       = qg_p        , ni        = ni_p        , nr         = nr_p         , &
                  pii      = pi_p        , p         = pres_p      , dz         = dz_p         , &
                  w        = w_p         , dt_in     = dt_microp   , itimestep  = itimestep    , &
                  rainnc   = rainnc_p    , rainncv   = rainncv_p   , snownc     = snownc_p     , &
                  snowncv  = snowncv_p   , graupelnc = graupelnc_p , graupelncv = graupelncv_p , &
                  sr       = sr_p        , rainprod  = rainprod_p  , evapprod   = evapprod_p   , &
                  re_cloud = recloud_p   , re_ice    = reice_p     , re_snow    = resnow_p     , &
                  has_reqc = has_reqc    , has_reqi  = has_reqi    , has_reqs   = has_reqs     , &
                  ntc      = ntc_p       , muc       = muc_p       , pflrain    = prcpfluxr_p  , &
                  pflice   = prcpfluxi_p , pflsnow   = prcpfluxs_p , pflgraul   = prcpfluxg_p  , &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde        , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme        , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte          &
                           )
       istep = istep + 1
       enddo
       call mpas_timer_stop('mp_thompson')

    case ("mp_thompson_aerosols")
       call mpas_timer_start('mp_thompson_aerosols')
       istep = 1
       do while (istep .le. n_microp)
       call mp_gt_driver( &
                  th       = th_p        , qv        = qv_p        , qc         = qc_p         , &
                  qr       = qr_p        , qi        = qi_p        , qs         = qs_p         , &
                  qg       = qg_p        , ni        = ni_p        , nr         = nr_p         , &
                  pii      = pi_p        , p         = pres_p      , dz         = dz_p         , &
                  w        = w_p         , dt_in     = dt_microp   , itimestep  = itimestep    , &
                  rainnc   = rainnc_p    , rainncv   = rainncv_p   , snownc     = snownc_p     , &
                  snowncv  = snowncv_p   , graupelnc = graupelnc_p , graupelncv = graupelncv_p , &
                  sr       = sr_p        , rainprod  = rainprod_p  , evapprod   = evapprod_p   , &
                  re_cloud = recloud_p   , re_ice    = reice_p     , re_snow    = resnow_p     , &
                  has_reqc = has_reqc    , has_reqi  = has_reqi    , has_reqs   = has_reqs     , &
                  nc       = nc_p        , nifa      = nifa_p      , nwfa       = nwfa_p       , &
                  nifa2d   = nifa2d_p    , nwfa2d    = nwfa2d_p    , ntc        = ntc_p        , &
                  muc      = muc_p       , pflrain   = prcpfluxr_p , pflice     = prcpfluxi_p  , &
                  pflsnow  = prcpfluxs_p , pflgraul  = prcpfluxg_p                             , &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde        , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme        , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte          &
                           )
          istep = istep + 1
          enddo
          call mpas_timer_stop('mp_thompson_aerosols')

    case ("mp_thompson_gocart2G")
       call mpas_timer_start('mp_thompson_gocart2G')
#ifdef DO_GOCART2G
! --- 1. Retrieve all the pointers from MPAS pools
       call mpas_pool_get_subpool(block_structs, 'NI2G_diags', NI2G_diags)
       call mpas_pool_get_subpool(block_structs, 'SU2G_diags', SU2G_diags)
       call mpas_pool_get_subpool(block_structs, 'CA2G_bc_diags', CA2G_bc_diags)
       call mpas_pool_get_subpool(block_structs, 'CA2G_br_diags', CA2G_br_diags)
       call mpas_pool_get_subpool(block_structs, 'CA2G_oc_diags', CA2G_oc_diags)
       call mpas_pool_get_subpool(block_structs, 'DU2G_diags', DU2G_diags)
       call mpas_pool_get_subpool(block_structs, 'SS2G_diags', SS2G_diags)

       call mpas_pool_get_array(NI2G_diags, 'niWT_bin1', niWT_bin1) 
       call mpas_pool_get_array(NI2G_diags, 'niWT_bin2', niWT_bin2)
       call mpas_pool_get_array(NI2G_diags, 'niWT_bin3', niWT_bin3) 
       call mpas_pool_get_array(NI2G_diags, 'nh3WT', nh3WT)
       call mpas_pool_get_array(NI2G_diags, 'nh4WT', nh4WT)
       call mpas_pool_get_array(SU2G_diags, 'suWT_so2', suWT_so2) 
       call mpas_pool_get_array(SU2G_diags, 'suWT_so4', suWT_so4)
       call mpas_pool_get_array(SU2G_diags, 'suWT_dms', suWT_dms) 
       call mpas_pool_get_array(SU2G_diags, 'suWT_msa', suWT_msa)
       call mpas_pool_get_array(CA2G_bc_diags, 'bcWT_phobic', bcWT_phobic) 
       call mpas_pool_get_array(CA2G_bc_diags, 'bcWT_philic', bcWT_philic)
       call mpas_pool_get_array(CA2G_br_diags, 'brWT_phobic', brWT_phobic)
       call mpas_pool_get_array(CA2G_br_diags, 'brWT_philic', brWT_philic)
       call mpas_pool_get_array(CA2G_oc_diags, 'ocWT_phobic', ocWT_phobic) 
       call mpas_pool_get_array(CA2G_oc_diags, 'ocWT_philic', ocWT_philic)
       call mpas_pool_get_array(DU2G_diags, 'duWT_bin1', duWT_bin1) 
       call mpas_pool_get_array(DU2G_diags, 'duWT_bin2', duWT_bin2)
       call mpas_pool_get_array(DU2G_diags, 'duWT_bin3', duWT_bin3) 
       call mpas_pool_get_array(DU2G_diags, 'duWT_bin4', duWT_bin4)
       call mpas_pool_get_array(DU2G_diags, 'duWT_bin5', duWT_bin5)
       call mpas_pool_get_array(SS2G_diags, 'ssWT_bin1', ssWT_bin1) 
       call mpas_pool_get_array(SS2G_diags, 'ssWT_bin2', ssWT_bin2)
       call mpas_pool_get_array(SS2G_diags, 'ssWT_bin3', ssWT_bin3) 
       call mpas_pool_get_array(SS2G_diags, 'ssWT_bin4', ssWT_bin4)
       call mpas_pool_get_array(SS2G_diags, 'ssWT_bin5', ssWT_bin5)

       offset = lbound(niWT_bin1, 1) - its

! --- 2. INITIALIZATION (Zero out garbage memory) ---
       if (itimestep == 1) then
          do i = its, ite
             idx = i + offset

             ! Nitrate and Ammonia
             niWT_bin1(idx) = 0.0_RKIND
             niWT_bin2(idx) = 0.0_RKIND
             niWT_bin3(idx) = 0.0_RKIND
             nh3WT(idx)     = 0.0_RKIND
             nh4WT(idx)     = 0.0_RKIND

             ! Sulfur Species
             suWT_so2(idx)  = 0.0_RKIND
             suWT_so4(idx)  = 0.0_RKIND
             suWT_dms(idx)  = 0.0_RKIND
             suWT_msa(idx)  = 0.0_RKIND

             ! Carbonaceous Aerosols
             bcWT_phobic(idx) = 0.0_RKIND
             bcWT_philic(idx) = 0.0_RKIND
             brWT_phobic(idx) = 0.0_RKIND
             brWT_philic(idx) = 0.0_RKIND
             ocWT_phobic(idx) = 0.0_RKIND
             ocWT_philic(idx) = 0.0_RKIND

             ! Mineral Dust
             duWT_bin1(idx) = 0.0_RKIND
             duWT_bin2(idx) = 0.0_RKIND
             duWT_bin3(idx) = 0.0_RKIND
             duWT_bin4(idx) = 0.0_RKIND
             duWT_bin5(idx) = 0.0_RKIND

             ! Sea Salt
             ssWT_bin1(idx) = 0.0_RKIND
             ssWT_bin2(idx) = 0.0_RKIND
             ssWT_bin3(idx) = 0.0_RKIND
             ssWT_bin4(idx) = 0.0_RKIND
             ssWT_bin5(idx) = 0.0_RKIND
          end do
       end if

! --- 3. LOAD THE STATE AND CALCULATE COLUMN MASS BEFORE (Add to Pointer) ---
       call mpas_chem_gocart2G%gocart2G_forMPASthom(state,time_lev)
          do ic = 1, mpas_chem_gocart2G%nchem
             do i = its, ite
                idx = i + offset
                col_mass_before = 0.0_RKIND
                do k = kts, kte
                   col_mass_before = col_mass_before + mpas_chem_gocart2G%chem_mr(i,k,ic) * rho_p(i,k,1) * dz_p(i,k,1)
                end do
                
                select case(ic)
                   case(1);     bcWT_phobic(idx) = bcWT_phobic(idx) + col_mass_before
                   case(2);     bcWT_philic(idx) = bcWT_philic(idx) + col_mass_before
                   case(3);     brWT_phobic(idx) = brWT_phobic(idx) + col_mass_before
                   case(4);     brWT_philic(idx) = brWT_philic(idx) + col_mass_before
                   case(5);     ocWT_phobic(idx) = ocWT_phobic(idx) + col_mass_before
                   case(6);     ocWT_philic(idx) = ocWT_philic(idx) + col_mass_before
                   case(7);     duWT_bin1(idx)   = duWT_bin1(idx)   + col_mass_before
                   case(8);     duWT_bin2(idx)   = duWT_bin2(idx)   + col_mass_before
                   case(9);     duWT_bin3(idx)   = duWT_bin3(idx)   + col_mass_before
                   case(10);    duWT_bin4(idx)   = duWT_bin4(idx)   + col_mass_before
                   case(11);    duWT_bin5(idx)   = duWT_bin5(idx)   + col_mass_before
                   case(12);    niWT_bin1(idx)   = niWT_bin1(idx)   + col_mass_before
                   case(13);    niWT_bin2(idx)   = niWT_bin2(idx)   + col_mass_before
                   case(14);    niWT_bin3(idx)   = niWT_bin3(idx)   + col_mass_before
                   case(15,16); suWT_so2(idx)    = suWT_so2(idx)    + col_mass_before
                   case(17,18); suWT_so4(idx)    = suWT_so4(idx)    + col_mass_before
                   case(19);    ssWT_bin1(idx)   = ssWT_bin1(idx)   + col_mass_before
                   case(20);    ssWT_bin2(idx)   = ssWT_bin2(idx)   + col_mass_before
                   case(21);    ssWT_bin3(idx)   = ssWT_bin3(idx)   + col_mass_before
                   case(22);    ssWT_bin4(idx)   = ssWT_bin4(idx)   + col_mass_before
                   case(23);    ssWT_bin5(idx)   = ssWT_bin5(idx)   + col_mass_before
                   case(24);    suWT_dms(idx)    = suWT_dms(idx)    + col_mass_before
                   case(25);    suWT_msa(idx)    = suWT_msa(idx)    + col_mass_before
                   case(26);    nh3WT(idx)       = nh3WT(idx)       + col_mass_before
                   case(27);    nh4WT(idx)       = nh4WT(idx)       + col_mass_before
                end select
             end do
          end do
 
! --- 4. EXECUTE THOMPSON MICROPHYSICS SEQUENCE ---
       do j = jts,jte
          do i = its,ite
             do k = kts,kte
                nwfa_p(i,k,j) = mpas_chem_gocart2G%chem_nwfa(i,k)
                nifa_p(i,k,j) = mpas_chem_gocart2G%chem_nifa(i,k)
             enddo
          enddo
       enddo
       istep = 1
       do while (istep .le. n_microp)
       call mp_gt_driver( &
                  th       = th_p        , qv        = qv_p        , qc         = qc_p         , &
                  qr       = qr_p        , qi        = qi_p        , qs         = qs_p         , &
                  qg       = qg_p        , ni        = ni_p        , nr         = nr_p         , &
                  pii      = pi_p        , p         = pres_p      , dz         = dz_p         , &
                  w        = w_p         , dt_in     = dt_microp   , itimestep  = itimestep    , &
                  rainnc   = rainnc_p    , rainncv   = rainncv_p   , snownc     = snownc_p     , &
                  snowncv  = snowncv_p   , graupelnc = graupelnc_p , graupelncv = graupelncv_p , &
                  sr       = sr_p        , rainprod  = rainprod_p  , evapprod   = evapprod_p   , &
                  re_cloud = recloud_p   , re_ice    = reice_p     , re_snow    = resnow_p     , &
                  has_reqc = has_reqc    , has_reqi  = has_reqi    , has_reqs   = has_reqs     , &
                  nc       = nc_p        , ntc       = ntc_p       , muc        = muc_p        , &
                  pflrain  = prcpfluxr_p , pflice    = prcpfluxi_p , pflsnow    = prcpfluxs_p  , &
                  pflgraul = prcpfluxg_p , nwfa      = nwfa_p      , nifa       = nifa_p       , &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde        , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme        , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte          &
                           )
       istep = istep + 1
       enddo
       do j = jts,jte
          call gs_scavenging_thermo(its, ite, kts, kte, &
                                    mpas_chem_gocart2G%nchem, dt_microp, itimestep, &
                                    th_p(:,:,j), pi_p(:,:,j), rho_p(:,:,j), & 
                                    dz_p(:,:,j), & ! <--- New: Passing layer thickness [m]
                                    qc_p(:,:,j), qi_p(:,:,j), &
                                    rainprod_p(:,:,j), rainncv_p(:,j), &
                                    mpas_chem_gocart2G%chem_mr, &
                                    mpas_chem_gocart2G%is_aerosol, &
                                    mpas_chem_gocart2G%scav_data_gas)
          do i = its,ite
             do k = kts,kte
                mpas_chem_gocart2G%chem_nwfa(i,k) = nwfa_p(i,k,j)
                mpas_chem_gocart2G%chem_nifa(i,k) = nifa_p(i,k,j)
             enddo
          enddo
       enddo
       call mpas_chem_gocart2G%gocart2G_fromMPASthom(state,time_lev)

! --- 5. CALCULATE COLUMN MASS AFTER (Subtract from Pointer) ---
          do ic = 1, mpas_chem_gocart2G%nchem
             do i = its, ite
                idx = i + offset
                col_mass_after = 0.0_RKIND
                do k = kts, kte
                   col_mass_after = col_mass_after + mpas_chem_gocart2G%chem_mr(i,k,ic) * rho_p(i,k,1) * dz_p(i,k,1)
                end do
                
                select case(ic)
                   case(1);     bcWT_phobic(idx) = bcWT_phobic(idx) - col_mass_after
                   case(2);     bcWT_philic(idx) = bcWT_philic(idx) - col_mass_after
                   case(3);     brWT_phobic(idx) = brWT_phobic(idx) - col_mass_after
                   case(4);     brWT_philic(idx) = brWT_philic(idx) - col_mass_after
                   case(5);     ocWT_phobic(idx) = ocWT_phobic(idx) - col_mass_after
                   case(6);     ocWT_philic(idx) = ocWT_philic(idx) - col_mass_after
                   case(7);     duWT_bin1(idx)   = duWT_bin1(idx)   - col_mass_after
                   case(8);     duWT_bin2(idx)   = duWT_bin2(idx)   - col_mass_after
                   case(9);     duWT_bin3(idx)   = duWT_bin3(idx)   - col_mass_after
                   case(10);    duWT_bin4(idx)   = duWT_bin4(idx)   - col_mass_after
                   case(11);    duWT_bin5(idx)   = duWT_bin5(idx)   - col_mass_after
                   case(12);    niWT_bin1(idx)   = niWT_bin1(idx)   - col_mass_after
                   case(13);    niWT_bin2(idx)   = niWT_bin2(idx)   - col_mass_after
                   case(14);    niWT_bin3(idx)   = niWT_bin3(idx)   - col_mass_after
                   case(15,16); suWT_so2(idx)    = suWT_so2(idx)    - col_mass_after
                   case(17,18); suWT_so4(idx)    = suWT_so4(idx)    - col_mass_after
                   case(19);    ssWT_bin1(idx)   = ssWT_bin1(idx)   - col_mass_after
                   case(20);    ssWT_bin2(idx)   = ssWT_bin2(idx)   - col_mass_after
                   case(21);    ssWT_bin3(idx)   = ssWT_bin3(idx)   - col_mass_after
                   case(22);    ssWT_bin4(idx)   = ssWT_bin4(idx)   - col_mass_after
                   case(23);    ssWT_bin5(idx)   = ssWT_bin5(idx)   - col_mass_after
                   case(24);    suWT_dms(idx)    = suWT_dms(idx)    - col_mass_after
                   case(25);    suWT_msa(idx)    = suWT_msa(idx)    - col_mass_after
                   case(26);    nh3WT(idx)       = nh3WT(idx)       - col_mass_after
                   case(27);    nh4WT(idx)       = nh4WT(idx)       - col_mass_after
                end select
             end do
          end do
#endif
       call mpas_timer_stop('mp_thompson_aerosols')

    case ("mp_wsm6")
       call mpas_timer_start('mp_wsm6')
       call wsm6( &
                  th        = th_p        , q          = qv_p         , qc        = qc_p      , &
                  qr        = qr_p        , qi         = qi_p         , qs        = qs_p      , &
                  qg        = qg_p        , den        = rho_p        , pii       = pi_p      , &
                  p         = pres_p      , delz       = dz_p         , delt      = dt_microp , &
                  g         = gravity     , cpd        = cp           , cpv       = cpv       , &
                  rd        = R_d         , rv         = R_v          , t0c       = svpt0     , &
                  ep1       = ep_1        , ep2        = ep_2         , qmin      = epsilon   , &
                  xls       = xls         , xlv0       = xlv          , xlf0      = xlf       , &
                  den0      = rho_a       , denr       = rho_w        , cliq      = cliq      , &
                  cice      = cice        , psat       = psat         , rain      = rainnc_p  , &
                  rainncv   = rainncv_p   , snow       = snownc_p     , snowncv   = snowncv_p , &
                  graupel   = graupelnc_p , graupelncv = graupelncv_p , sr        = sr_p      , &
                  re_cloud  = recloud_p   , re_ice     = reice_p      , re_snow   = resnow_p  , &
                  has_reqc  = has_reqc    , has_reqi   = has_reqi     , has_reqs  = has_reqs  , &
                  re_qc_bg  = re_qc_bg    , re_qi_bg   = re_qi_bg     , re_qs_bg  = re_qs_bg  , &
                  re_qc_max = re_qc_max   , re_qi_max  = re_qi_max    , re_qs_max = re_qs_max , &
                  errmsg    = errmsg      , errflg     = errflg                               , &
                  ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde       , &
                  ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme       , &
                  its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte         &
               )
       call mpas_timer_stop('mp_wsm6')

    case default
 end select microp_select

!... calculate the 10cm radar reflectivity and relative humidity, if needed:
 if (l_diags) then
    !ensure that we only call compute_radar_reflectivity() if we are using an MPS that supports
    !the computation of simulated radar reflectivity:
    if(trim(microp_scheme) == "mp_wsm6"              .or. &
       trim(microp_scheme) == "mp_thompson"          .or. &
       trim(microp_scheme) == "mp_thompson_aerosols" .or. &
       trim(microp_scheme) == "mp_thompson_gocart2G") then
       call compute_radar_reflectivity(configs,diag_physics,its,ite)
    else
       call mpas_log_write('*** NOTICE: NOT computing simulated radar reflectivity')
       call mpas_log_write('            since WSM6 or Thompson microphysics scheme was not selected')
    end if

    !calculate the relative humidity over water if the temperature is strictly greater than 0.C,
    !over ice otherwise.
    call compute_relhum(diag,its,ite)
 end if

!... copy updated precipitation from the wrf-physics grid back to the geodesic-dynamics grid:
 call precip_to_MPAS(configs,diag_physics,its,ite)

!... copy updated cloud microphysics variables from the wrf-physics grid back to the geodesic-
!    dynamics grid:
 call microphysics_to_MPAS(configs,mesh,state,time_lev,diag,diag_physics,tend_physics,tend,its,ite)

!... deallocation of all microphysics arrays:
!$OMP BARRIER
!$OMP MASTER
 call deallocate_microphysics(configs)
!$OMP END MASTER

!call mpas_log_write('---enter subroutine driver_microphysics:')
!call mpas_log_write('')

 end subroutine driver_microphysics

!=================================================================================================================
 subroutine precip_from_MPAS(configs,diag_physics,its,ite)
!=================================================================================================================

!input variables:
 type(mpas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!output variables:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 integer,pointer:: nCellsSolve
 real,dimension(:),pointer:: graupelncv,rainncv,snowncv,sr

!local variables and arrays:
 integer:: i,j

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

 call mpas_pool_get_array(diag_physics,'graupelncv',graupelncv)
 call mpas_pool_get_array(diag_physics,'rainncv'   ,rainncv   )
 call mpas_pool_get_array(diag_physics,'snowncv'   ,snowncv   )
 call mpas_pool_get_array(diag_physics,'sr'        ,sr        )

!variables common to all cloud microphysics schemes:
 do j = jts, jte
 do i = its, ite
    rainncv_p(i,j) = 0._RKIND
    rainnc_p(i,j)  = 0._RKIND
 enddo
 enddo

 do i = its,ite
    rainncv(i) = 0._RKIND
 enddo

!variables specific to different cloud microphysics schemes:
 microp_select: select case(trim(microp_scheme))
    case ("mp_thompson","mp_thompson_aerosols","mp_thompson_gocart2G","mp_wsm6")
       do j = jts, jte
       do i = its, ite
          snowncv_p(i,j)    = 0._RKIND
          graupelncv_p(i,j) = 0._RKIND
          snownc_p(i,j)     = 0._RKIND
          graupelnc_p(i,j)  = 0._RKIND
          sr_p(i,j)         = 0._RKIND
       enddo
       enddo

       do i = its,ite
          snowncv(i)    = 0._RKIND
          graupelncv(i) = 0._RKIND
          sr(i)         = 0._RKIND
       enddo

    case default
 end select microp_select

 end subroutine precip_from_MPAS

!=================================================================================================================
 subroutine precip_to_MPAS(configs,diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 integer,dimension(:),pointer:: i_rainnc

 real(kind=RKIND),pointer:: config_bucket_rainnc
 real(kind=RKIND),dimension(:),pointer:: precipw
 real(kind=RKIND),dimension(:),pointer:: graupelnc,rainnc,snownc
 real(kind=RKIND),dimension(:),pointer:: graupelncv,rainncv,snowncv,sr

!local variables and arrays:
 integer:: i,j,k
 real(kind=RKIND):: rho_a

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme       )
 call mpas_pool_get_config(configs,'config_bucket_rainnc',config_bucket_rainnc)

 call mpas_pool_get_array(diag_physics,'i_rainnc'  ,i_rainnc  )
 call mpas_pool_get_array(diag_physics,'precipw'   ,precipw   )
 call mpas_pool_get_array(diag_physics,'graupelnc' ,graupelnc )
 call mpas_pool_get_array(diag_physics,'graupelncv',graupelncv)
 call mpas_pool_get_array(diag_physics,'rainnc'    ,rainnc    )
 call mpas_pool_get_array(diag_physics,'rainncv'   ,rainncv   )
 call mpas_pool_get_array(diag_physics,'snownc'    ,snownc    )
 call mpas_pool_get_array(diag_physics,'snowncv'   ,snowncv   )
 call mpas_pool_get_array(diag_physics,'sr'        ,sr        )

 do i = its,ite
    precipw(i) = 0._RKIND
 enddo

!variables common to all cloud microphysics schemes:
 do j = jts,jte
 do i = its,ite

    !precipitable water:
    do k = kts,kte
       rho_a = rho_p(i,k,j) / (1._RKIND + qv_p(i,k,j))
       precipw(i) = precipw(i) + qv_p(i,k,j) * rho_a * dz_p(i,k,j)
    enddo

    !time-step precipitation:
    rainncv(i) = rainnc_p(i,j)

    !accumulated precipitation:
    rainnc(i) = rainnc(i) + rainncv(i)

   if(l_acrain .and. config_bucket_rainnc.gt.0._RKIND .and. &
      rainnc(i).gt.config_bucket_rainnc) then
      i_rainnc(i) = i_rainnc(i) + 1
      rainnc(i)   = rainnc(i) - config_bucket_rainnc
   endif

 enddo
 enddo

!variables specific to different cloud microphysics schemes:
 microp_select: select case(trim(microp_scheme))
    case ("mp_thompson","mp_thompson_aerosols","mp_thompson_gocart2G","mp_wsm6")
       do j = jts,jte
       do i = its,ite
          !time-step precipitation:
          snowncv(i)    = snownc_p(i,j)
          graupelncv(i) = graupelnc_p(i,j)
          sr(i) = (snownc_p(i,j) + graupelnc_p(i,j)) / (rainnc_p(i,j)+1.e-12)

          !accumulated precipitation:
          snownc(i)    = snownc(i) + snowncv(i)
          graupelnc(i) = graupelnc(i) + graupelncv(i)
       enddo
       enddo

    case default
 end select microp_select

 end subroutine precip_to_MPAS

!=================================================================================================================
 subroutine compute_radar_reflectivity(configs,diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 character(len=StrKIND),pointer:: microp_scheme
 real(kind=RKIND),dimension(:,:),pointer:: refl10cm
 real(kind=RKIND),dimension(:),pointer:: refl10cm_max,refl10cm_1km,refl10cm_1km_max

!local variables and arrays:
 integer:: i,j,k,kp
 real(kind=RKIND),dimension(:),allocatable:: qv1d,qc1d,qr1d,qs1d,qg1d,t1d,p1d,nr1d,dBZ1d,zp
 real(kind=RKIND):: w1,w2

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme',microp_scheme)

 call mpas_pool_get_array(diag_physics,'refl10cm',refl10cm)
 call mpas_pool_get_array(diag_physics,'refl10cm_max',refl10cm_max)
 call mpas_pool_get_array(diag_physics,'refl10cm_1km',refl10cm_1km)
 call mpas_pool_get_array(diag_physics,'refl10cm_1km_max',refl10cm_1km_max)

 microp_select: select case(trim(microp_scheme))
    case ("mp_kessler")
       call physics_error_fatal('--- calculation of radar reflectivity is not available' // &
                                 'with kessler cloud microphysics')

    case ("mp_wsm6")
       if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
       if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
       if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
       if(.not.allocated(qr1d) ) allocate(qr1d(kts:kte) )
       if(.not.allocated(qs1d) ) allocate(qs1d(kts:kte) )
       if(.not.allocated(qg1d) ) allocate(qg1d(kts:kte) )
       if(.not.allocated(dBz1d)) allocate(dBZ1d(kts:kte))
       if(.not.allocated(zp)   ) allocate(zp(kts:kte)   )

       do j = jts,jte
       do i = its,ite
          do k = kts,kte
             p1d(k) = pres_p(i,k,j)
             t1d(k) = th_p(i,k,j) * pi_p(i,k,j)
             qv1d(k)  = qv_p(i,k,j)
             qr1d(k)  = qr_p(i,k,j)
             qs1d(k)  = qs_p(i,k,j)
             qg1d(k)  = qg_p(i,k,j)
             dBZ1d(k) = -35._RKIND
             zp(k) = z_p(i,k,j) - z_p(i,1,j) + 0.5*dz_p(i,k,j) ! height AGL
          enddo

          call refl10cm_wsm6(qv1d,qr1d,qs1d,qg1d,t1d,p1d,dBZ1d,kts,kte)

          kp = 1
          do k = kts,kte
             dBZ1d(k) = max(-35._RKIND,dBZ1d(k))
             refl10cm(k,i) = dBZ1d(k)
             if(zp(k) .lt. 1000.) kp = k
          enddo
          refl10cm_max(i) = maxval(dBZ1d(:))
          w1 = (zp(kp+1)-1000.)/(zp(kp+1)-zp(kp))
          w2 = 1.0 - w1
          refl10cm_1km(i) = w1*dBZ1d(kp) + w2*dBZ1d(kp+1)
          refl10cm_1km_max(i) = max(refl10cm_1km_max(i),refl10cm_1km(i))
       enddo
       enddo

       if(allocated(p1d)  ) deallocate(p1d  )
       if(allocated(t1d)  ) deallocate(t1d  )
       if(allocated(qv1d) ) deallocate(qv1d )
       if(allocated(qr1d) ) deallocate(qr1d )
       if(allocated(qs1d) ) deallocate(qs1d )
       if(allocated(qg1d) ) deallocate(qg1d )
       if(allocated(dBz1d)) deallocate(dBZ1d)
       if(allocated(zp)   ) deallocate(zp   )

    case ("mp_thompson","mp_thompson_aerosols","mp_thompson_gocart2G")
       if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
       if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
       if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
       if(.not.allocated(qc1d) ) allocate(qc1d(kts:kte) )
       if(.not.allocated(qr1d) ) allocate(qr1d(kts:kte) )
       if(.not.allocated(qs1d) ) allocate(qs1d(kts:kte) )
       if(.not.allocated(qg1d) ) allocate(qg1d(kts:kte) )
       if(.not.allocated(nr1d) ) allocate(nr1d(kts:kte) )
       if(.not.allocated(dBz1d)) allocate(dBZ1d(kts:kte))
       if(.not.allocated(zp)   ) allocate(zp(kts:kte)   )

       do j = jts,jte
       do i = its,ite
          do k = kts,kte
             p1d(k) = pres_p(i,k,j)
             t1d(k) = th_p(i,k,j) * pi_p(i,k,j)
             qv1d(k)  = qv_p(i,k,j)
             qc1d(k)  = qc_p(i,k,j)
             qr1d(k)  = qr_p(i,k,j)
             qs1d(k)  = qs_p(i,k,j)
             qg1d(k)  = qg_p(i,k,j)
             nr1d(k)  = nr_p(i,k,j)
             dBZ1d(k) = -35._RKIND
             zp(k) = z_p(i,k,j) - z_p(i,1,j) + 0.5*dz_p(i,k,j) ! height AGL
          enddo

          call calc_refl10cm(qv1d,qc1d,qr1d,nr1d,qs1d,qg1d,t1d,p1d,dBZ1d,kts,kte,i,j)

          kp = 1
          do k = kts,kte
             dBZ1d(k) = max(-35._RKIND,dBZ1d(k))
             refl10cm(k,i) = dBZ1d(k)
             if(zp(k) .lt. 1000.) kp = k
          enddo
          refl10cm_max(i) = maxval(dBZ1d(:))
          w1 = (zp(kp+1)-1000.)/(zp(kp+1)-zp(kp))
          w2 = 1.0 - w1
          refl10cm_1km(i) = w1*dBZ1d(kp) + w2*dBZ1d(kp+1)
          refl10cm_1km_max(i) = max(refl10cm_1km_max(i),refl10cm_1km(i))
       enddo
       enddo

       if(allocated(p1d)  ) deallocate(p1d  )
       if(allocated(t1d)  ) deallocate(t1d  )
       if(allocated(qv1d) ) deallocate(qv1d )
       if(allocated(qc1d) ) deallocate(qc1d )
       if(allocated(qr1d) ) deallocate(qr1d )
       if(allocated(qs1d) ) deallocate(qs1d )
       if(allocated(qg1d) ) deallocate(qg1d )
       if(allocated(nr1d) ) deallocate(nr1d )
       if(allocated(dBz1d)) deallocate(dBZ1d)
       if(allocated(zp)   ) deallocate(zp   )

    case default
 end select microp_select

 end subroutine compute_radar_reflectivity

!=================================================================================================================
 subroutine compute_relhum(diag,its,ite)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag

!local pointers:
 real(kind=RKIND),dimension(:,:),pointer:: relhum

!local variables and arrays:
 integer:: i,j,k

 real(kind=RKIND):: tempc, rh
 real(kind=RKIND),dimension(:),allocatable:: qv1d,qvs1d,t1d,p1d

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_array(diag,'relhum',relhum)
 do k = kts,kte
    do i = its,ite
       relhum(k,i) = 0._RKIND
    enddo
 enddo

 if(.not.allocated(p1d)  ) allocate(p1d(kts:kte)  )
 if(.not.allocated(t1d)  ) allocate(t1d(kts:kte)  )
 if(.not.allocated(qv1d) ) allocate(qv1d(kts:kte) )
 if(.not.allocated(qvs1d)) allocate(qvs1d(kts:kte))

 do j = jts,jte
 do i = its,ite
    do k = kts,kte
       p1d(k)  = pres_p(i,k,j)
       t1d(k)  = th_p(i,k,j) * pi_p(i,k,j)
       tempc   = t1d(k) - 273.16_RKIND
       qvs1d(k) = rslf(p1d(k),t1d(k))
       if(tempc .le. 0._RKIND) qvs1d(k) = rsif(p1d(k),t1d(k))
       qv1d(k) = qv_p(i,k,j)
       relhum(k,i) = qv1d(k) / qvs1d(k)
       relhum(k,i) = relhum(k,i) * 100._RKIND
    enddo
 enddo
 enddo

 if(allocated(p1d)  ) deallocate(p1d  )
 if(allocated(t1d)  ) deallocate(t1d  )
 if(allocated(qv1d) ) deallocate(qv1d )
 if(allocated(qvs1d)) deallocate(qvs1d)

 end subroutine compute_relhum

!=================================================================================================================
     subroutine gs_scavenging_thermo(its, ite, kts, kte, nchem, dt, itimestep, &
                                     th_p, pi_p, rho, dz, qc, qi, rainprod, rainncv, &
                                     chem, is_aerosol, scav_data_gas)
!=================================================================================================================
! Comprehensive Scavenging Routine
! 1. Hydrophobic Aerosols: Laakso et al. (2003) empirical washout
! 2. Gases: Dynamic T and dz dependent Henry's Law scavenging
!=================================================================================================================
        implicit none
        integer, intent(in) :: its, ite, kts, kte, nchem, itimestep
        real(kind=RKIND), intent(in) :: dt
        real(kind=RKIND), dimension(its:ite, kts:kte), intent(in) :: th_p, pi_p, rho, dz, qc, qi, rainprod
        real(kind=RKIND), dimension(its:ite), intent(in) :: rainncv
        real(kind=RKIND), dimension(its:ite, kts:kte, nchem), intent(inout) :: chem
        logical, dimension(nchem), intent(in) :: is_aerosol
        real(kind=RKIND), dimension(6, nchem), intent(in) :: scav_data_gas

        ! Local variables
        integer :: i, k, jn
        real(kind=RKIND) :: t_local, tfac, h_phys, h_star, k1, k2, h_ion, scav_eff
        real(kind=RKIND) :: lambda, lambda_in, lambda_out, conv_factor
        real(kind=RKIND) :: rain_mm_hr, total_water_vol
        
        ! Constants
        real(kind=RKIND), parameter :: r_gas_const = 0.08206_RKIND  ! L*atm/(mol*K)
        real(kind=RKIND), parameter :: laakso_a    = 1.0e-4_RKIND   ! s-1
        real(kind=RKIND), parameter :: laakso_b    = 0.6_RKIND      ! unitless
        real(kind=RKIND), parameter :: t_ref       = 298.15_RKIND   ! Kelvin

        logical :: found_activity = .false. 

        h_ion = 1.0e-5_RKIND ! pH 5.0 (Standard for cloud/rainwater)

        do i = its, ite
            ! rainncv [mm] / dt [s] * 3600 [s/hr] = [mm/hr]
            rain_mm_hr = (rainncv(i) / dt) * 3600.0_RKIND

            do k = kts, kte
                t_local = th_p(i,k) * pi_p(i,k)
                
                ! Convert mixing ratio (kg/kg) to volumetric mass (kg/m3)
                ! total_water_vol units: kg/m3
                total_water_vol = max(1.0e-7_RKIND, (qc(i,k) + qi(i,k)) * rho(i,k))

                do jn = 1, nchem

                    ! ---------------------------------------------------------
                    ! PHASE A: AEROSOL SCAVENGING
                    ! ---------------------------------------------------------
                    if (is_aerosol(jn)) then
                        ! Apply Laakso (2003) to Hydrophobic portions: BC(1), OC(3), BrC(5)
                        if (jn == 1 .or. jn == 3 .or. jn == 5) then
                            if (rain_mm_hr > 0.01_RKIND) then
                                lambda = laakso_a * (rain_mm_hr**laakso_b)
                            else
                                lambda = 0.0_RKIND
                            end if
                        else
                            cycle ! Thompson handles hydrophilic/activated aerosols
                        end if

                    ! ---------------------------------------------------------
                    ! PHASE B: GAS-PHASE SCAVENGING (Dynamic Thermo)
                    ! ---------------------------------------------------------
                    else if (scav_data_gas(1, jn) > 0.0_RKIND) then
                        
                        ! 1. Thermodynamics (Van't Hoff Temperature Correction)
                        tfac = (1.0_RKIND / t_local) - (1.0_RKIND / t_ref)
                        h_phys = scav_data_gas(1, jn) * exp(scav_data_gas(2, jn) * tfac)
                        
                        ! Calculate Effective Henry's Law (H*) for dissociation
                        if (scav_data_gas(3, jn) > 0.0_RKIND) then
                            k1 = scav_data_gas(3, jn) * exp(scav_data_gas(4, jn) * tfac)
                            if (scav_data_gas(5, jn) > 0.0_RKIND) then
                                k2 = scav_data_gas(5, jn) * exp(scav_data_gas(6, jn) * tfac)
                                if (scav_data_gas(1, jn) > 10.0_RKIND) then
                                    h_star = h_phys * (1.0_RKIND + (k1 * h_ion) / k2) ! NH3
                                else
                                    h_star = h_phys * (1.0_RKIND + k1/h_ion + (k1*k2)/(h_ion**2)) ! SO2
                                end if
                            else
                                h_star = h_phys * (1.0_RKIND + k1/h_ion) ! MSA
                            end if
                        else
                            h_star = h_phys ! DMS / Insoluble
                        end if

                        ! 2. Log-Linear Solubility Efficiency for Rainout (In-cloud)
                        if (h_star <= 0.1_RKIND) then
                            scav_eff = 0.0_RKIND
                        else if (h_star >= 1.0e6_RKIND) then
                            scav_eff = 1.0_RKIND
                        else
                            scav_eff = (log10(max(1.0e-10_RKIND, h_star)) + 1.0_RKIND) / 7.0_RKIND
                        end if

                        ! 3. Calculate Total Scavenging Rate [1/s]
                        ! lambda_in: Consumption of gas by in-situ precipitation production
                        lambda_in = (rainprod(i,k) / total_water_vol) * scav_eff
                        
                        ! lambda_out: Dynamic washout based on RT/dz equilibrium
                        ! conversion factor = (R_gas * T) / (3600 * 1000)
                        conv_factor = (r_gas_const * t_local) / 3600000.0_RKIND
                        lambda_out = rain_mm_hr * h_star * (conv_factor / dz(i,k))

                        lambda = lambda_in + lambda_out

                    else
                        cycle ! No data for this gas
                    end if

                    ! ---------------------------------------------------------
                    ! PHASE C: APPLY MASS REMOVAL
                    ! ---------------------------------------------------------
                    if (lambda > 0.0_RKIND) then
                        chem(i,k,jn) = chem(i,k,jn) * exp(-min(lambda * dt, 10.0_RKIND))

                        ! Diagnostic Output for first signs of activity
                         !if (itimestep <= 50 .and. (.not. found_activity) .and. lambda > 1.e-7_RKIND) then
                         !   call mpas_log_write('GS_SCAV_ACTUAL [Step $i, Col $i, Lev $i, Spec $i]', &
                         !       intArgs=(/itimestep, i, k, jn/))
                         !   call mpas_log_write('   -> RainRate=$r, H*=$r, DZ=$r, Lambda=$r', &
                         !       realArgs=(/rain_mm_hr, h_star, dz(i,k), lambda/))
                         !   found_activity = .true. 
                         !end if 
                    end if

                end do ! jn loop
            end do ! k loop
        end do ! i loop
      end subroutine gs_scavenging_thermo

!=================================================================================================================
 end module mpas_atmphys_driver_microphysics
!=================================================================================================================
