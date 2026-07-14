!=================================================================================================================
 module module_cu_ntiedtke
 use mpas_kind_types,only: kind_phys => RKIND
 use mpas_log
 use cu_ntiedtke,only: cu_ntiedtke_run, &
                       cu_ntiedtke_init
 use cu_ntiedtke_common
 use cu_ntiedtke_post,only: cu_ntiedtke_post_run
 use cu_ntiedtke_pre,only: cu_ntiedtke_pre_run

#ifdef DO_GOCART2G
 use mpas_chemistry_gocart2G_vars
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_dimension, mpas_pool_get_subpool
 use cu_ntiedtke_ctrans_post
 use cu_ntiedtke_ctrans_pre
 use cu_ntiedtke_ctrans
#endif


 implicit none
 private
 public:: cu_ntiedtke_driver, &
          ntiedtkeinit


 contains


!=================================================================================================================
 subroutine cu_ntiedtke_driver(                          &
                 dt,itimestep,stepcu                     &
                ,raincv,pratec,qfx,hfx                   &
                ,u3d,v3d,w,t3d,qv3d,qc3d,qi3d,pi3d,rho3d &
                ,qvften,thften                           &
                ,dz8w,pcps,p8w,xland,cu_act_flag,dx      &
                ,f_qc,f_qi                               &
                ,grav,xlf,xls,xlv,rd,rv,cp               &
                ,rthcuten,rqvcuten,rqccuten,rqicuten     &
                ,rucuten,rvcuten                         &
                ,do_ctrchems,do_scavchems                &
                ,block_structs                           &
                ,ids,ide,jds,jde,kds,kde                 &
                ,ims,ime,jms,jme,kms,kme                 &
                ,its,ite,jts,jte,kts,kte                 &
                ,errmsg,errflg)
!=================================================================================================================
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- th3d        3d potential temperature (k)
!-- t3d         temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- qc3d        3d cloud mixing ratio (kg/kg)
!-- qi3d        3d ice mixing ratio (kg/kg)
!-- rho3d       3d air density (kg/m^3)
!-- p8w         3d hydrostatic pressure at full levels (pa)
!-- pcps        3d hydrostatic pressure at half levels (pa)
!-- pi3d        3d exner function (dimensionless)
!-- qvften      3d total advective + PBL moisture tendency (kg kg-1 s-1)
!-- thften      3d total advective + PBL + radiative temperature tendency (k s-1)
!-- rthcuten      theta tendency due to 
!                 cumulus scheme precipitation (k/s)
!-- rucuten       u wind tendency due to 
!                 cumulus scheme precipitation (k/s)
!-- rvcuten       v wind tendency due to 
!                 cumulus scheme precipitation (k/s)
!-- rqvcuten      qv tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- rqccuten      qc tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- rqicuten      qi tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- rainc         accumulated total cumulus scheme precipitation (mm)
!-- raincv        cumulus scheme precipitation (mm)
!-- pratec        precipitiation rate from cumulus scheme (mm/s)
!-- dz8w        dz between full levels (m)
!-- qfx         upward moisture flux at the surface (kg/m^2/s)
!-- hfx         upward heat flux at the surface (w/m^2) 
!-- dt          time step (s)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-----------------------------------------------------------------------------------------------------------------

!--- input arguments:
 logical,intent(in):: do_ctrchems,do_scavchems
 logical,intent(in),optional:: f_qc,f_qi

 integer,intent(in):: ids,ide,jds,jde,kds,kde, &
                      ims,ime,jms,jme,kms,kme, &
                      its,ite,jts,jte,kts,kte

 integer,intent(in):: itimestep,stepcu

 real(kind=kind_phys),intent(in):: cp,grav,rd,rv,xlf,xls,xlv

 real(kind=kind_phys),intent(in):: dt

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: dx,hfx,qfx,xland

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    dz8w,   &
    pcps,   &
    p8w,    &
    pi3d,   &
    qc3d,   &
    qvften, &
    thften, &
    qi3d,   &
    qv3d,   &
    rho3d,  &
    t3d,    &
    u3d,    &
    v3d,    &
    w

!--- inout arguments:
 logical,intent(inout),dimension(ims:ime,jms:jme):: cu_act_flag

 real(kind=kind_phys),intent(inout),dimension(ims:ime,jms:jme):: raincv, pratec

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme),optional:: &
 rqccuten,  &
 rqicuten,  &
 rqvcuten,  &
 rthcuten,  &
 rucuten,   &
 rvcuten

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg

!--- local variables and arrays:
 integer:: i,im,j,k,kx,kx1
 integer,dimension(its:ite)::  slimsk

 real(kind=kind_phys):: delt
 real(kind=kind_phys),dimension(its:ite):: rn
 real(kind=kind_phys),dimension(its:ite,kts:kte):: prsl,omg,ghtl
 real(kind=kind_phys),dimension(its:ite,kts:kte):: uf,vf,tf,qvf,qcf,qif
 real(kind=kind_phys),dimension(its:ite,kts:kte):: qvftenz,thftenz
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: prsi,ghti,zi

 real(kind=kind_phys),dimension(its:ite):: dx_hv,hfx_hv,qfx_hv,xland_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: dz_hv,pi_hv,prsl_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: qv_hv,qc_hv,qi_hv,rho_hv,t_hv,u_hv,v_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: qvften_hv,thften_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: prsi_hv,w_hv

 real(kind=kind_phys),dimension(its:ite):: raincv_hv,pratec_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: rthcuten_hv,rqvcuten_hv,rqccuten_hv,rqicuten_hv, &
                                                   rucuten_hv,rvcuten_hv

!--- additional local arrays output from the subroutine cu_ntiedtke_run:
!    aerosol species:
 logical,dimension(its:ite):: ldcum,lddraf
 integer,dimension(its:ite):: kctype,kcbot,kctop,kdtop
 real(kind=kind_phys),dimension(its:ite,kts:kte):: pmfu,pmfd
 real(kind=kind_phys),dimension(its:ite,kts:kte):: pmfu_entr,pmfd_entr,pmfu_detr,pmfd_detr

 integer:: ic,nchem
 real(kind=kind_phys),dimension(:),allocatable:: fscav
 real(kind=kind_phys),dimension(:,:,:),allocatable:: chem,chemf
 real(kind=kind_phys),dimension(:,:,:),allocatable:: ptenc,ptencf

 type(mpas_pool_type), intent(inout) :: block_structs
 type (mpas_pool_type), pointer :: CA2G_bc_diags, CA2G_br_diags, CA2G_oc_diags
 type (mpas_pool_type), pointer :: DU2G_diags, NI2G_diags, SS2G_diags, SU2G_diags 

 logical, dimension(:), allocatable:: is_aerosol, is_coarse
 real(kind=kind_phys), dimension(:,:), allocatable:: rate_incloud, rate_washout
 real(kind=kind_phys), dimension(:), pointer :: bcSV_phobic,bcSV_philic 
 real(kind=kind_phys), dimension(:), pointer :: brSV_phobic,brSV_philic
 real(kind=kind_phys), dimension(:), pointer :: ocSV_phobic,ocSV_philic
 real(kind=kind_phys), dimension(:), pointer :: duSV_bin1,duSV_bin2,duSV_bin3,duSV_bin4,duSV_bin5
 real(kind=kind_phys), dimension(:), pointer :: niSV_bin1, niSV_bin2, niSV_bin3,nh3SV,nh4SV
 real(kind=kind_phys), dimension(:), pointer :: ssSV_bin1,ssSV_bin2,ssSV_bin3,ssSV_bin4,ssSV_bin5
 real(kind=kind_phys), dimension(:), pointer :: suSV_dms,suSV_so2,suSV_so4,suSV_msa

! --- Manual Indexing Sequence from gocart2G_forMPASphys_init ---
 integer, parameter :: n_bc_phobic = 1,  n_bc_philic = 2
 integer, parameter :: n_br_phobic = 3,  n_br_philic = 4
 integer, parameter :: n_oc_phobic = 5,  n_oc_philic = 6
 integer, parameter :: n_du1 = 7, n_du2 = 8, n_du3 = 9, n_du4 = 10, n_du5 = 11
 integer, parameter :: n_ni1 = 12, n_ni2 = 13, n_ni3 = 14
 integer, parameter :: n_so2 = 15, n_so2_vol = 16
 integer, parameter :: n_so4 = 17, n_so4_vol = 18
 integer, parameter :: n_ss1 = 19, n_ss2 = 20, n_ss3 = 21, n_ss4 = 22, n_ss5 = 23
 integer, parameter :: n_dms = 24, n_msa = 25
 integer, parameter :: n_nh3 = 26, n_nh4 = 27
 integer :: idx, offset 

!-----------------------------------------------------------------------------------------------------------------

 errmsg = ' '
 errflg = 0

#ifdef DO_GOCART2G
! Retrieve all the diagnostic pools
     call mpas_pool_get_subpool(block_structs, 'CA2G_bc_diags', CA2G_bc_diags)
     call mpas_pool_get_subpool(block_structs, 'CA2G_br_diags', CA2G_br_diags)
     call mpas_pool_get_subpool(block_structs, 'CA2G_oc_diags', CA2G_oc_diags)
     call mpas_pool_get_subpool(block_structs, 'DU2G_diags',    DU2G_diags)
     call mpas_pool_get_subpool(block_structs, 'NI2G_diags',    NI2G_diags)
     call mpas_pool_get_subpool(block_structs, 'SS2G_diags',    SS2G_diags)
     call mpas_pool_get_subpool(block_structs, 'SU2G_diags',    SU2G_diags)

     call mpas_pool_get_array(CA2G_bc_diags,'bcSV_phobic',bcSV_phobic)
     call mpas_pool_get_array(CA2G_bc_diags,'bcSV_philic',bcSV_philic)
    
     call mpas_pool_get_array(CA2G_br_diags,'brSV_phobic',brSV_phobic)
     call mpas_pool_get_array(CA2G_br_diags,'brSV_philic',brSV_philic)
    
     call mpas_pool_get_array(CA2G_oc_diags,'ocSV_phobic',ocSV_phobic)
     call mpas_pool_get_array(CA2G_oc_diags,'ocSV_philic',ocSV_philic)
    
     call mpas_pool_get_array(DU2G_diags,'duSV_bin1',duSV_bin1)
     call mpas_pool_get_array(DU2G_diags,'duSV_bin2',duSV_bin2)
     call mpas_pool_get_array(DU2G_diags,'duSV_bin3',duSV_bin3)
     call mpas_pool_get_array(DU2G_diags,'duSV_bin4',duSV_bin4)
     call mpas_pool_get_array(DU2G_diags,'duSV_bin5',duSV_bin5)
    
     call mpas_pool_get_array(NI2G_diags,'niSV_bin1',niSV_bin1)
     call mpas_pool_get_array(NI2G_diags,'niSV_bin2',niSV_bin2)
     call mpas_pool_get_array(NI2G_diags,'niSV_bin3',niSV_bin3)
     call mpas_pool_get_array(NI2G_diags,'nh3SV',nh3SV)
     call mpas_pool_get_array(NI2G_diags,'nh4SV',nh4SV)
    
     call mpas_pool_get_array(SS2G_diags,'ssSV_bin1',ssSV_bin1)
     call mpas_pool_get_array(SS2G_diags,'ssSV_bin2',ssSV_bin2)
     call mpas_pool_get_array(SS2G_diags,'ssSV_bin3',ssSV_bin3)
     call mpas_pool_get_array(SS2G_diags,'ssSV_bin4',ssSV_bin4)
     call mpas_pool_get_array(SS2G_diags,'ssSV_bin5',ssSV_bin5)
    
     call mpas_pool_get_array(SU2G_diags,'suSV_dms',suSV_dms)
     call mpas_pool_get_array(SU2G_diags,'suSV_so2',suSV_so2)
     call mpas_pool_get_array(SU2G_diags,'suSV_so4',suSV_so4)
     call mpas_pool_get_array(SU2G_diags,'suSV_msa',suSV_msa)

! --- Initializing Integrated Sinks to Zero at Start of Simulation ---
     offset = lbound(bcSV_phobic, 1) - its 
     if (itimestep == 1) then
        do i = its, ite
           idx = i + offset       
           ! Carbonaceous
           bcSV_phobic(idx) = 0.0_kind_phys
           bcSV_philic(idx) = 0.0_kind_phys
           brSV_phobic(idx) = 0.0_kind_phys
           brSV_philic(idx) = 0.0_kind_phys
           ocSV_phobic(idx) = 0.0_kind_phys
           ocSV_philic(idx) = 0.0_kind_phys
 
           ! Dust bins
           duSV_bin1(idx)   = 0.0_kind_phys
           duSV_bin2(idx)   = 0.0_kind_phys
           duSV_bin3(idx)   = 0.0_kind_phys
           duSV_bin4(idx)   = 0.0_kind_phys
           duSV_bin5(idx)   = 0.0_kind_phys
 
           ! Nitrate and Ammonia
           niSV_bin1(idx)   = 0.0_kind_phys
           niSV_bin2(idx)   = 0.0_kind_phys
           niSV_bin3(idx)   = 0.0_kind_phys
           nh3SV(idx)       = 0.0_kind_phys
           nh4SV(idx)       = 0.0_kind_phys
 
           ! Sea salt bins
           ssSV_bin1(idx)   = 0.0_kind_phys
           ssSV_bin2(idx)   = 0.0_kind_phys
           ssSV_bin3(idx)   = 0.0_kind_phys
           ssSV_bin4(idx)   = 0.0_kind_phys
           ssSV_bin5(idx)   = 0.0_kind_phys
 
           ! Sulfur group
           suSV_dms(idx)    = 0.0_kind_phys
           suSV_so2(idx)    = 0.0_kind_phys
           suSV_so4(idx)    = 0.0_kind_phys
           suSV_msa(idx)    = 0.0_kind_phys
        end do
     end if
#endif

 call cu_ntiedtke_init( &
                       con_cp  = cp  , con_rd = rd   , con_rv   = rv   , con_xlv = xlv    , &
                       con_xls = xls , con_xlf = xlf , con_grav = grav , errmsg  = errmsg , &
                       errflg  = errflg                                                     &
                      )

 do j = jts,jte
    do i = its,ite
       cu_act_flag(i,j)=.true.
    enddo
 enddo

 do j = jts,jte

    do i = its,ite
       dx_hv(i)    = dx(i,j)
       hfx_hv(i)   = hfx(i,j)
       qfx_hv(i)   = qfx(i,j)
       xland_hv(i) = xland(i,j)
    enddo

    do k = kts,kte
       do i = its,ite
          dz_hv(i,k)   = dz8w(i,k,j)
          pi_hv(i,k)   = pi3d(i,k,j)
          prsl_hv(i,k) = pcps(i,k,j)
          qv_hv(i,k)   = qv3d(i,k,j)
          qc_hv(i,k)   = qc3d(i,k,j)
          qi_hv(i,k)   = qi3d(i,k,j)
          rho_hv(i,k)  = rho3d(i,k,j)
          t_hv(i,k)    = t3d(i,k,j)
          u_hv(i,k)    = u3d(i,k,j)
          v_hv(i,k)    = v3d(i,k,j)

          qvften_hv(i,k) = qvften(i,k,j)
          thften_hv(i,k) = thften(i,k,j)
       enddo
    enddo
    do k = kts,kte+1
       do i = its,ite
          prsi_hv(i,k) = p8w(i,k,j)
          w_hv(i,k)    = w(i,k,j)
       enddo
    enddo

    call cu_ntiedtke_pre_run( &
                     its    = its       , ite     = ite     , kts     = kts     , kte       = kte       , &
                     im     = im        , kx      = kx      , kx1     = kx1     , itimestep = itimestep , &
                     stepcu = stepcu    , dt      = dt      , grav    = grav    , xland     = xland_hv  , &
                     dz     = dz_hv     , pres    = prsl_hv , presi   = prsi_hv , t         = t_hv      , &
                     rho    = rho_hv    , qv      = qv_hv   , qc      = qc_hv   , qi        = qi_hv     , &
                     u      = u_hv      , v       = v_hv    , w       = w_hv    , qvften    = qvften_hv , &
                     thften = thften_hv , qvftenz = qvftenz , thftenz = thftenz , slimsk    = slimsk    , &
                     delt   = delt      , prsl    = prsl    , ghtl    = ghtl    , tf        = tf        , &
                     qvf    = qvf       , qcf     = qcf     , qif     = qif     , uf        = uf        , &
                     vf     = vf        , prsi    = prsi    , ghti    = ghti    , omg       = omg       , &
                     errmsg = errmsg    , errflg  = errflg                                                &
                            )

    call cu_ntiedtke_run( &
                     pu        =  uf       , pv        = vf        , pt        = tf        , &
                     pqv       = qvf       , pqc       = qcf       , pqi       = qif       , &
                     pqvf      = qvftenz   , ptf       = thftenz   , poz       = ghtl      , &
                     pzz       = ghti      , pomg      = omg       , pap       = prsl      , &
                     paph      = prsi      , evap      = qfx_hv    , hfx       = hfx_hv    , &
                     zprecc    = rn        , lndj      = slimsk    , lq        = im        , &
                     km        = kx        , km1       = kx1       , dt        = delt      , &
                     dx        = dx_hv     , ldcum     = ldcum     , lddraf    = lddraf    , &
                     kctype    = kctype    , kcbot     = kcbot     , kctop     = kctop     , &
                     kdtop     = kdtop     , pmfu      = pmfu      , pmfd      = pmfd      , &
                     pmfu_entr = pmfu_entr , pmfd_entr = pmfd_entr , pmfu_detr = pmfu_detr , &
                     pmfd_detr = pmfd_detr , errmsg    = errmsg    , errflg    = errflg      &
                        )

    call cu_ntiedtke_post_run( &
                     its       = its         , ite       = ite         , kts       = kts         , &
                     kte       = kte         , stepcu    = stepcu      , dt        = dt          , &
                     exner     = pi_hv       , qv        = qv_hv       , qc        = qc_hv       , &
                     qi        = qi_hv       , t         = t_hv        , u         = u_hv        , &
                     v         = v_hv        , qvf       = qvf         , qcf       = qcf         , &
                     qif       = qif         , tf        = tf          , uf        = uf          , &
                     vf        = vf          , rn        = rn          , raincv    = raincv_hv   , &
                     pratec    = pratec_hv   , rthcuten  = rthcuten_hv , rqvcuten  = rqvcuten_hv , &
                     rqccuten  = rqccuten_hv , rqicuten  = rqicuten_hv , rucuten   = rucuten_hv  , &
                     rvcuten   = rvcuten_hv  , errmsg    = errmsg      , errflg    = errflg        &
                             )

    do i = its,ite
       raincv(i,j) = raincv_hv(i)
       pratec(i,j) = pratec_hv(i)
    enddo

    do k = kts,kte
       do i = its,ite
          rucuten(i,k,j)  = rucuten_hv(i,k)
          rvcuten(i,k,j)  = rvcuten_hv(i,k)
          rthcuten(i,k,j) = rthcuten_hv(i,k)
          rqvcuten(i,k,j) = rqvcuten_hv(i,k)
       enddo
    enddo

    if(present(rqccuten))then
       if(f_qc) then
          do k = kts,kte
             do i = its,ite
                rqccuten(i,k,j) = rqccuten_hv(i,k)
             enddo
          enddo
       endif
    endif

    if(present(rqicuten))then
       if(f_qi) then
          do k = kts,kte
             do i = its,ite
                rqicuten(i,k,j) = rqicuten_hv(i,k)
             enddo
          enddo
       endif
    endif

#ifdef DO_GOCART2G
    nchem = mpas_chem_gocart2G%nchem
    if(do_ctrchems) then
       if(.not.allocated(fscav) ) allocate(fscav(nchem))
       if(.not.allocated(is_aerosol)) allocate(is_aerosol(nchem))
       if(.not.allocated(is_coarse)) allocate(is_coarse(nchem))
       if(.not.allocated(rate_incloud)) allocate(rate_incloud(its:ite,nchem))
       if(.not.allocated(rate_washout)) allocate(rate_washout(its:ite,nchem))
       if(.not.allocated(chem)  ) allocate(chem(its:ite,kts:kte,nchem)  )
       if(.not.allocated(chemf) ) allocate(chemf(its:ite,kts:kte,nchem) )
       if(.not.allocated(ptenc) ) allocate(ptenc(its:ite,kts:kte,nchem) )
       if(.not.allocated(ptencf)) allocate(ptencf(its:ite,kts:kte,nchem))
! --- Initialize Tile Array for this step ---
       rate_incloud(:,:) = 0.0_kind_phys
       rate_washout(:,:) = 0.0_kind_phys

       do ic = 1,nchem
          fscav(ic) = mpas_chem_gocart2G%fscav(ic)
          chem(:,:,ic) = mpas_chem_gocart2G%chem_mr(:,:,ic)
          is_aerosol(ic) = mpas_chem_gocart2G%is_aerosol(ic)
          is_coarse(ic) = mpas_chem_gocart2G%is_coarse(ic)
       enddo
       call cu_ntiedtke_ctrans_pre_run( &
               its    = its    , ite  = ite  , kts   = kts   , kte    = kte    , &
               nchem  = nchem  , chem = chem , chemf = chemf , errmsg = errmsg , &
               errflg = errflg                                                   &
                                      )

       call cu_ntiedtke_ctrans_run( &
               klon          = im                 , &
               klev          = kx                 , &
               nchem         = nchem              , &
               ldcum         = ldcum              , &
               lddraf        = lddraf             , &
               kctype        = kctype             , &
               kcbot         = kcbot              , &
               kctop         = kctop              , &
               kdtop         = kdtop              , &
               grav          = grav               , &
               ztmst         = delt               , &
               do_scav       = do_scavchems       , &
               fscav         = fscav              , &
               is_aerosol    = mpas_chem_gocart2G%is_aerosol, &
               scav_data_gas = mpas_chem_gocart2G%scav_data_gas, &
               chem          = chemf              , &
               thv           = t_hv               , &
               ptenc         = ptenc              , &
               ghti          = ghti               , &
               paph          = prsi               , &
               pmfu          = pmfu               , &
               pmfd          = pmfd               , &
               pmfude_rate   = pmfu_detr          , &
               pmfdde_rate   = pmfd_detr          , &
               errmsg        = errmsg             , &
               errflg        = errflg             , &
               itimestep     = itimestep          , &
               rate_incloud  = rate_incloud         &
                                   )

       call cu_ntiedtke_ctrans_post_run( &
               its    = its   , ite   = ite   , kts    = kts    , kte    = kte    , &
               nchem  = nchem , ptenc = ptenc , ptencf = ptencf , errmsg = errmsg , &
               errflg = errflg &
                                       )
       if (do_scavchems) then
          call cu_ntiedtke_below_cloud_scavenge(its, ite, kts, kte, nchem, delt, itimestep, &
                                                kcbot, pratec_hv, prsi_hv, t_hv, dz_hv, rho_hv, chem, is_aerosol, &
                                                is_coarse, mpas_chem_gocart2G%scav_data_gas, rate_washout, ptencf)
       endif

! --- Final Transfer to Diagnostic Pointers (Cumulative kg/m2) ---
     do i = its, ite
      if (i == its) then
         call mpas_log_write('DEBUG SCAV [Step $i]: delt=$r, PTR=$r, RateIn=$r, RateWash=$r', &
                        intArgs=(/ itimestep /), &
                        realArgs=(/ delt, bcSV_phobic(i), &
                                   rate_incloud(i, n_bc_phobic), &
                                   rate_washout(i, n_bc_phobic) /))
      end if
        idx = i + offset
        ! Black Carbon
        bcSV_phobic(idx) = bcSV_phobic(idx) + ((rate_incloud(i, n_bc_phobic) + rate_washout(i, n_bc_phobic)) * delt)
        bcSV_philic(idx) = bcSV_philic(idx) + ((rate_incloud(i, n_bc_philic) + rate_washout(i, n_bc_philic)) * delt)
    
        ! Brown Carbon
        brSV_phobic(idx) = brSV_phobic(idx) + ((rate_incloud(i, n_br_phobic) + rate_washout(i, n_br_phobic)) * delt)
        brSV_philic(idx) = brSV_philic(idx) + ((rate_incloud(i, n_br_philic) + rate_washout(i, n_br_philic)) * delt)
    
        ! Organic Carbon
        ocSV_phobic(idx) = ocSV_phobic(idx) + ((rate_incloud(i, n_oc_phobic) + rate_washout(i, n_oc_phobic)) * delt)
        ocSV_philic(idx) = ocSV_philic(idx) + ((rate_incloud(i, n_oc_philic) + rate_washout(i, n_oc_philic)) * delt)
    
        ! Mineral Dust
        duSV_bin1(idx) = duSV_bin1(idx) + ((rate_incloud(i, n_du1) + rate_washout(i, n_du1)) * delt)
        duSV_bin2(idx) = duSV_bin2(idx) + ((rate_incloud(i, n_du2) + rate_washout(i, n_du2)) * delt)
        duSV_bin3(idx) = duSV_bin3(idx) + ((rate_incloud(i, n_du3) + rate_washout(i, n_du3)) * delt)
        duSV_bin4(idx) = duSV_bin4(idx) + ((rate_incloud(i, n_du4) + rate_washout(i, n_du4)) * delt)
        duSV_bin5(idx) = duSV_bin5(idx) + ((rate_incloud(i, n_du5) + rate_washout(i, n_du5)) * delt)
    
        ! Nitrate Bins
        niSV_bin1(idx) = niSV_bin1(idx) + ((rate_incloud(i, n_ni1) + rate_washout(i, n_ni1)) * delt)
        niSV_bin2(idx) = niSV_bin2(idx) + ((rate_incloud(i, n_ni2) + rate_washout(i, n_ni2)) * delt)
        niSV_bin3(idx) = niSV_bin3(idx) + ((rate_incloud(i, n_ni3) + rate_washout(i, n_ni3)) * delt)
    
        ! Sulfur Species (Summing volcanic and non-volcanic into one diagnostic)
        suSV_so2(idx) = suSV_so2(idx) + ((rate_incloud(i, n_so2) + rate_washout(i, n_so2)) * delt) + &
                                   ((rate_incloud(i, n_so2_vol) + rate_washout(i, n_so2_vol)) * delt)
        
        suSV_so4(idx) = suSV_so4(idx) + ((rate_incloud(i, n_so4) + rate_washout(i, n_so4)) * delt) + &
                                   ((rate_incloud(i, n_so4_vol) + rate_washout(i, n_so4_vol)) * delt)
        
        suSV_dms(idx) = suSV_dms(idx) + ((rate_incloud(i, n_dms) + rate_washout(i, n_dms)) * delt)
        suSV_msa(idx) = suSV_msa(idx) + ((rate_incloud(i, n_msa) + rate_washout(i, n_msa)) * delt)
    
        ! Sea Salt Bins
        ssSV_bin1(idx) = ssSV_bin1(idx) + ((rate_incloud(i, n_ss1) + rate_washout(i, n_ss1)) * delt)
        ssSV_bin2(idx) = ssSV_bin2(idx) + ((rate_incloud(i, n_ss2) + rate_washout(i, n_ss2)) * delt)
        ssSV_bin3(idx) = ssSV_bin3(idx) + ((rate_incloud(i, n_ss3) + rate_washout(i, n_ss3)) * delt)
        ssSV_bin4(idx) = ssSV_bin4(idx) + ((rate_incloud(i, n_ss4) + rate_washout(i, n_ss4)) * delt)
        ssSV_bin5(idx) = ssSV_bin5(idx) + ((rate_incloud(i, n_ss5) + rate_washout(i, n_ss5)) * delt)
    
        ! Ammonia and Ammonium
        nh3SV(idx) = nh3SV(idx) + ((rate_incloud(i, n_nh3) + rate_washout(i, n_nh3)) * delt)
        nh4SV(idx) = nh4SV(idx) + ((rate_incloud(i, n_nh4) + rate_washout(i, n_nh4)) * delt)
     end do

       do ic = 1,nchem
          mpas_chem_gocart2G%chemcuten(:,:,ic) = ptencf(:,:,ic)
       enddo

       if(allocated(fscav) ) deallocate(fscav)
       if(allocated(chem)  ) deallocate(chem  )
       if(allocated(chemf) ) deallocate(chemf )
       if(allocated(ptenc) ) deallocate(ptenc )
       if(allocated(ptencf)) deallocate(ptencf)
       if(allocated(is_aerosol)) deallocate(is_aerosol)
       if(allocated(is_coarse)) deallocate(is_coarse)
       if(allocated(rate_incloud)) deallocate(rate_incloud)
       if(allocated(rate_washout)) deallocate(rate_washout)
    else
       do ic = 1,nchem
          mpas_chem_gocart2G%chemcuten(:,:,ic) = 0._kind_phys
       enddo
    endif
#endif

 enddo

 end subroutine cu_ntiedtke_driver

!=================================================================================================================
 subroutine ntiedtkeinit(rthcuten,rqvcuten,rqccuten,rqicuten, &
                         rucuten,rvcuten,rthften,rqvften,     &
                         restart,p_qc,p_qi,p_first_scalar,    &
                         allowed_to_read,                     &
                         ids, ide, jds, jde, kds, kde,        &
                         ims, ime, jms, jme, kms, kme,        &
                         its, ite, jts, jte, kts, kte)
!=================================================================================================================

!--- input arguments:
 logical,intent(in):: allowed_to_read,restart

 integer,intent(in):: ids, ide, jds, jde, kds, kde, &
                      ims, ime, jms, jme, kms, kme, &
                      its, ite, jts, jte, kts, kte
 integer,intent(in):: p_first_scalar,p_qi,p_qc

!--- output arguments:
 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme )::  &
    rthcuten,rqvcuten,rqccuten,rqicuten,rucuten,rvcuten,rthften,rqvften

!--- local variables and arrays:
 integer:: i,j,k,itf,jtf,ktf

!-----------------------------------------------------------------------------------------------------------------

 jtf = min0(jte,jde-1)
 ktf = min0(kte,kde-1)
 itf = min0(ite,ide-1)

 if(.not.restart)then
    do j = jts,jtf
       do k = kts,ktf
          do i = its,itf
             rthcuten(i,k,j) = 0.
             rqvcuten(i,k,j) = 0.
             rucuten(i,k,j)  = 0.
             rvcuten(i,k,j)  = 0.
          enddo
       enddo
    enddo

    do j = jts,jtf
       do k = kts,ktf
          do i = its,itf
             rthften(i,k,j)=0.
             rqvften(i,k,j)=0.
          enddo
       enddo
    enddo

    if(p_qc .ge. p_first_scalar) then
       do j = jts,jtf
          do k = kts,ktf
             do i = its,itf
                rqccuten(i,k,j)=0.
             enddo
          enddo
       enddo
    endif

    if(p_qi .ge. p_first_scalar) then
       do j = jts,jtf
          do k = kts,ktf
             do i = its,itf
                rqicuten(i,k,j)=0.
             enddo
          enddo
       enddo
    endif
 endif

 end subroutine ntiedtkeinit

!=================================================================================================================
subroutine cu_ntiedtke_below_cloud_scavenge(its, ite, kts, kte, nchem, dt, itimestep, &
                                            kcbot, pratec, prsi, t_hv, dz, rho, chem, is_aerosol, &
                                            is_coarse, scav_data_gas, rate_washout, ptencf)
!=================================================================================================================
! CONVECTIVE BELOW-CLOUD SCAVENGING (V8.1 PRODUCTION - METADATA DRIVEN)
! 1. Aerosols: Size-split scavenging (Laakso fine mode vs. UFS coarse mode via is_coarse flag)
! 2. Gases: Dynamic T and dz dependent Henry's Law scavenging
!=================================================================================================================
        implicit none

        ! Argument declarations
        integer, intent(in) :: its, ite, kts, kte, nchem, itimestep
        real(kind=kind_phys), intent(in) :: dt
        integer, dimension(its:ite), intent(in) :: kcbot
        real(kind=kind_phys), dimension(its:ite), intent(in) :: pratec
        real(kind=kind_phys), dimension(its:ite, kts:kte+1), intent(in) :: prsi
        real(kind=kind_phys), dimension(its:ite, kts:kte), intent(in) :: t_hv  ! Temperature [K]
        real(kind=kind_phys), dimension(its:ite, kts:kte), intent(in) :: dz, rho
        real(kind=kind_phys), dimension(its:ite, kts:kte, nchem), intent(inout) :: chem, ptencf
        real(kind=kind_phys), dimension(its:ite, nchem), intent(inout) :: rate_washout
        logical, dimension(nchem), intent(in) :: is_aerosol
        logical, dimension(nchem), intent(in) :: is_coarse ! <--- New: Clean size classification flag  

        real(kind=kind_phys), dimension(6, nchem), intent(in) :: scav_data_gas

        ! Local variables
        integer :: i, k, jn, k_bottomup_base 
        real(kind=kind_phys) :: rain_rate_mm_hr, lambda, chem_before
        real(kind=kind_phys) :: tfac, h_phys, h_star, k1, k2, h_ion, conv_factor

        ! Constants
        real(kind=kind_phys), parameter :: r_gas_const = 0.08206_kind_phys
        real(kind=kind_phys), parameter :: laakso_a    = 1.0e-4_kind_phys
        real(kind=kind_phys), parameter :: laakso_b    = 0.6_kind_phys
        real(kind=kind_phys), parameter :: t_ref       = 298.15_kind_phys

        h_ion = 1.0e-5_kind_phys
      
        do i = its, ite
            ! Convert surface convective rain rate to mm/hr
            rain_rate_mm_hr = pratec(i) * 3600.0_kind_phys
            
            if (rain_rate_mm_hr > 0.01_kind_phys) then
                
                ! kcbot is top-down; convert to bottom-up index for scavenging
                k_bottomup_base = kte - kcbot(i) + 1 

                do k = kts, k_bottomup_base
                    do jn = 1, nchem
                        lambda = 0.0_kind_phys

                        ! ---------------------------------------------------------
                        ! PHASE A: AEROSOL SCAVENGING (SIZE-SEGREGATED BRANCH)
                        ! ---------------------------------------------------------
                        if (is_aerosol(jn)) then

                            if (.not. is_coarse(jn)) then
                                ! --- FINE MODE BRANCH (Phase-Dependent) ---
                                if (t_hv(i,k) >= 268.15_kind_phys) then
                                    ! Liquid rain regime: Baseline Laakso (2003)
                                    lambda = laakso_a * (rain_rate_mm_hr**laakso_b)
                                else
                                    ! Ice/Snow regime: Amplified UFS-style structural scavenging
                                    ! (3.4 * 1.0e-4 = 3.4e-4) derived from coarse-mode coeffs (8.72e-4/2.56e-4 = 3.4)
                                    lambda = 3.4e-4_kind_phys * (rain_rate_mm_hr**0.96_kind_phys)
                                end if
                            else
                                ! --- COARSE MODE: UFS/GEOS-Chem Inertial Impaction Logic ---
                                if (t_hv(i,k) >= 268.15_kind_phys) then
                                    ! Liquid rain collection regime
                                    lambda = 2.56e-4_kind_phys * (rain_rate_mm_hr**0.79_kind_phys)
                                else
                                    ! Ice/Snow collection regime
                                    lambda = 8.72e-4_kind_phys * (rain_rate_mm_hr**0.96_kind_phys)
                                end if
                            end if

                        ! ---------------------------------------------------------
                        ! PHASE B: GAS SCAVENGING (DYNAMIC THERMO)
                        ! ---------------------------------------------------------
                        else if (scav_data_gas(1, jn) > 0.0_kind_phys) then
                            
                            tfac = (1.0_kind_phys / t_hv(i,k)) - (1.0_kind_phys / t_ref)
                            h_phys = scav_data_gas(1, jn) * exp(scav_data_gas(2, jn) * tfac)
                            
                            if (scav_data_gas(3, jn) > 0.0_kind_phys) then
                                k1 = scav_data_gas(3, jn) * exp(scav_data_gas(4, jn) * tfac)
                                if (scav_data_gas(5, jn) > 0.0_kind_phys) then
                                    k2 = scav_data_gas(5, jn) * exp(scav_data_gas(6, jn) * tfac)
                                    if (scav_data_gas(1, jn) > 10.0_kind_phys) then
                                        h_star = h_phys * (1.0_kind_phys + (k1 * h_ion) / k2) ! NH3 logic
                                    else
                                        h_star = h_phys * (1.0_kind_phys + k1/h_ion + (k1*k2)/(h_ion**2)) ! Acid logic
                                    end if
                                else
                                    h_star = h_phys * (1.0_kind_phys + k1/h_ion) ! Single dissociation
                                end if
                            else
                                h_star = h_phys ! Insoluble / No dissociation
                            end if

                            ! Calculate dynamic conversion factor (R*T / 3.6e6) and apply local dz
                            conv_factor = (r_gas_const * t_hv(i,k)) / 3600000.0_kind_phys
                            lambda = rain_rate_mm_hr * h_star * (conv_factor / dz(i,k))

                        end if
                        
                        ! ---------------------------------------------------------
                        ! PHASE C: APPLY MASS REMOVAL
                        ! ---------------------------------------------------------
                        if (lambda > 0.0_kind_phys) then
                            chem_before = chem(i,k,jn)
                            
                            ! Exponential decay with numerical limit
                            chem(i,k,jn) = chem(i,k,jn) * exp(-min(lambda * dt, 10.0_kind_phys))
                            
                            ! Update physics tendency for MPAS
                            ptencf(i,k,jn) = ptencf(i,k,jn) + (chem(i,k,jn) - chem_before) / dt
                            
                            ! Update surface diagnostic (kg/m2/s)
                            rate_washout(i, jn) = rate_washout(i, jn) + &
                                                  ((chem_before - chem(i,k,jn)) * rho(i,k) * dz(i,k)) / dt
                        end if
                    end do
                end do
            end if
        end do

     end subroutine cu_ntiedtke_below_cloud_scavenge

!=================================================================================================================
 end module module_cu_ntiedtke
!=================================================================================================================
