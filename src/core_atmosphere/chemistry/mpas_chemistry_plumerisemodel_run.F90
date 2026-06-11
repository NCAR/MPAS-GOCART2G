! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
! Run Plume Rise Model to Calculate Smoke Plume Injection Height
! Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-04-20.
!
! For reference: 
! https://github.com/wrf-model/WRF/blob/master/chem/module_chem_plumerise_scalar.F
! https://github.com/saulorfreitas/brams_SRF/blob/master/src/brams/ccatt/chem_plumerise_scalar.f90
!
!-------------------------------------------------------------------------------------------!
! Plume rise model for vegetation fires (CPTEC/INPE 2005-2006,2009)			            !
! Refs.:										                                     !
! Freitas, S. R., K. M. Longo, J. Trentmann, D. Latham. Technical Note: Sensitivity         !
! of 1D smoke plume rise models to the inclusion of environmental wind drag.                !
! Atmospheric Chemistry  and Physics, 2010.                                                 !
!                                                                                           !
! Freitas, S. R., K. M. Longo, R. Chatfield, D. Latham, M. A. F. Silva Dias, M. O. Andreae, !
! E. Prins, J. C. Santos, R. Gielow and J. A. Carvalho Jr.: Including the sub-grid scale    !
! plume rise of vegetation fires in low resolution atmospheric transport models. 	       !
!  Atmospheric Chemistry and Physics,2007.				                                !
! 											                                     !
! Freitas, S. R.; Longo, K. M.; M. Andreae. Impact of including the plume rise of vegetation! 
! fires in numerical simulations of associated atmospheric pollutants. Geophys. Res. Lett., !
! 33, L17808, doi:10.1029/2006GL026608, 2006.                                               !
!                                                                                           !
! Questions/bugs/comments:                                                                  !
! Saulo Freitas (saulo.freitas@inpe.br, saulo.r.de.freitas@gmail.com) 	                 !
! Jaqueline Pereira (jaqueline.pereira@inpe.br, pereira.njaqueline@gmail.com)               !
!-------------------------------------------------------------------------------------------!

!==================================================================================================================
 module mpas_chemistry_plumerisemodel_run

 use mpas_log
 use mpas_derived_types 
 use mpas_kind_types,only: kind_phys => RKIND

 implicit none
    
     private                      
     public:: plumerisemodel_run


     !------------------------------------------------------------------------------------------------------------
     ! Global configurations

          logical :: do_windeffect = .true.       !- turn ON (true), OFF (false) the env wind effect on plume rise
          integer :: sound_type = 0               !- set it to zero for the coupling with 3d host model. set 1 for reading a precreated sound already defined at PRM vertical grid
          !integer, parameter :: heat_flag = 1    !- set 1 for use prescribed heat flux values. set 2 to use heat flux from FRP
          logical :: use_zv_mtt = .true.          !- set true for MTT solution. set false for turn ON the virtual boundary height from Tarshish and Romps (2022)
          logical :: use_bott_adv = .true.        !- set true for using Bott's adv scheme. set false to use original adv scheme
          integer :: lflag = 1                    !- if using Bott's scheme: set the desired order (0-8).0 for 1st order. 1 for 2nd order, 2 for 4a. order...
          logical :: do_microphysics = .true.     !- turn ON, OFF cloud microphysics parameterization


     ! Global physical constants

          real(kind=kind_phys), parameter :: &
               grav  = 9.80665                  , & ![m/s2] av. gravitational acceleration
               cpd   = 1004.64                  , & ! [J/K/kg] specific heat at constant pressure
               rd    = 287.04                   , & !> [J/K/kg] gas constant
               pi    = 3.14159265358979323846264338327950288
        
     ! Global variables and parameters    
     
          integer, parameter :: nveg_ag =   4   !upper index of vegetation categories
          integer, parameter ::   nkp   = 200   !upper index of PRM vertical levels
          integer, parameter ::   ntime = 200   !
        
          real(kind=kind_phys), parameter :: alpha_mtt = 0.1   !- entrainment constant - range 0.1 - 0.16 Carazzo et al 2006


          real(kind=kind_phys), dimension(nkp) ::  &
               w,      & !
               t,      & !
               theq,   & !
               qv,     & !
               qc,     & !
               qh,     & !
               qi,     & !
               sc,     & !  
               vth,    & !
               vti,    & !
               rho,    & !
               txs,    & !
               est,    & !
               dqsdz,  & !
               qsat,   & !
               qpas,   & !
               qtotal, & !
               td,     & !
               vel_p,  & !
               rad_p,  & !
               alpha
          real(kind=kind_phys), dimension(nkp) ::  &
               wc,    & !
               wt,    & ! 
               tt,    & !
               qvt,   & !
               qct,   & !
               qht,   & !
               qit,   & !
               sct,   & !
               vel_t, & !
               rad_t    !
          real(kind=kind_phys), dimension(nkp) ::  &
               dzm,    & ! used as inverse of full level vertical grid (inverse of interface of the layer)
               dzt,    & ! used as inverse of half level vertical grid (inverse of midle of the layer)
               vctr1,  & !
               vctr2,  & !
               vt3dc,  & !
               vt3df,  & !
               vt3dk,  & !
               vt3dg,  & !
               scr1,   & !
               vt3dj,  & !
               vt3dn,  & !
               vt3do,  & !
               vt3da,  & !
               scr2,   & ! 
               vt3db,  & !
               vt3dd,  & !
               vt3dl,  & !
               vt3dm,  & !
               vt3di,  & !
               vt3dh,  & !
               vt3de,  & !
               rbuoy,  & !
               dwdt_entr !
          real(kind=kind_phys), dimension(nkp) ::  & ! environment at plume grid
               the,  & !
               thee, & !
               rhe,  & !
               sce,  & !
               tde     ! 
          real(kind=kind_phys), dimension(nkp) :: &
               te,    & ! temperature interpolated at PRM grid (kelvin)
               pe,    & ! pressure interpolated at PRM grid (KPa)
               upe,   & ! zonal wind interpolated at PRM grid (m.s-1)
               vpe,   & ! meridinal wind interpolated at PRM grid (m.s-1)
               qvenv, & ! environment water vapor mixing ratio interpolated at PRM grid (kg h20v/kg dry air)
               vel_e, & ! environmental horizontal wind computed from u and v at PRM grid(m.s-1)
               zt,    & ! thermo and water levels in the PRM grid (midle of the layer) (m above surface)
               zm       ! dynamical levels in the PRM model (interface of the layer) (m above surface)
          real(kind=kind_phys), dimension(1200) :: & ! environment at RAMS  grid
               wcon,   & !
               thtcon, & !
               zzcon,  & !  1D cartesian height of interface of layer of the host model (meters above surface)
               scon      !
          real(kind=kind_phys) :: &
               dz,        & !
               visc(nkp), & !
               viscosity, & !
               tstpf        !   
          integer :: &
               n,   & !
               nm1, & !
               l      !counter
          real(kind=kind_phys) :: & ! advection
               advw,      & !
               advt,      & !
               advv,      & !
               advc,      & !
               advh,      & !
               advi,      & !
               cvh(nkp),  & !
               cvi(nkp),  & !
               adiabat,   & !
               wbar,      & !
               vhrel,     & !
               virel        !  
          real(kind=kind_phys) :: &
               zsurf, & ! set to 0 in zero_plumegen_coms and never changed afterwards
               zbase, & !
               ztop     !
          integer :: LBASE
          real(kind=kind_phys) :: & ! entrain
               area,     & !
               rsurf,    & ! entrainment surface radius (m) set in get_fire_properties
               radius(nkp) !  
          real(kind=kind_phys) :: & ! heating
               heating(ntime), & ! heating rates computed in get_fire_properties
               fmoist,         & ! fuel moisture fraction set to 1/10 in get_fire_properties
               bload,          & ! total loading set to 10 in get_fire_properties
               heat_fluxW        ! cache variable for heatflux in get_fire_properties [W.m-2]
          real(kind=kind_phys) :: &
               dt,   & !
               time, & !
               tdur    ! number of seconds in the burn set in get_fire_properties to 60*mdur
          integer :: &
               mintime, & ! elapsed time, minutes - makeplume
               mdur,    & ! duration of burn, minutes 
               maxtime    ! set in in get_fire_properties
          real(kind=kind_phys) :: ztop_(ntime) !
    
          integer :: &
               kmt,              & ! index for the PRM model upper vertical grid. will be set to nkp in the subroutine get_env_condition
               imm,              & ! counter for lower and upper boundary of the heat flux
               k1,               & ! index for the lowest injection height when imm=1
               k2                  ! index for the highest injection height when imm=2

          real(kind=kind_phys) :: zbot_plume, ztop_plume !

     !------------------------------------------------------------------------------------------------------------


 contains

!=================================================================================================================
 subroutine plumerisemodel_run( &

               !variables from MPAS-A:
               its,ite, jts,jte, kts,kte                , &
               temp_in                                  , &
               pres_in                                  , &
               u_in                                     , &
               v_in                                     , &
               hml_in                                   , &
               qv_in                                    , &
               z_in                                     , &

               !variables to PRM (input)
               iveg_ag_in                               , & 
               firesize_avg_corr                        , &
               firesize_std_corr                        , &
               frp_avg_corr                             , &
               frp_std_corr                             , &
               heat_flag                                , &

               !variables from PRM (output)
               k1_zbot_out                              , &
               k2_ztop_out                              , &
               zbot_out                                 , &
               ztop_out                                   &                                
               )


     !*******************************************************************************
     implicit none
     

          !---------------------------------------------------------------
          ! Variables and arrays from/to plumerismodel_driver 
          

               integer,  intent(in   )   ::     &
                    its,ite,                    & ! index for the host model horizontal grid (its=1, ite = nCellsSolve)
                    jts,jte,                    & ! index for local horizontal grid (jts=jte=1)
                    kts,kte                       ! index for the host model vertical grid (kts=1, kte = nVertLevels)

               real(kind=kind_phys), dimension( its:ite, kts:kte, jts:jte ) ::     &
                    temp_in,                    & ! 3D Temperature of the host model (k)interpolated to theta points
                    pres_in,                    & ! 3D Pressure of the host model (Pa) interpolated to theta points
                    u_in,                       & ! 3D zonal wind of the host model (m/s) interpolated to theta points 
                    v_in,                       & ! 3D meridional wind of the host model (m/s) interpolated to theta points 
                    hml_in,                     & ! 3D geometric height of middle of layer of the host model (m MSL = meters above mean sea level)
                    qv_in,                      & ! 3D water vapor mixing ratio of the host model (kg h20vapor/kg dry air)
                    z_in                          ! 3D geometric height of layer interfaces of the host model located at w points (m MSL = meters above mean sea level)

               
               integer, dimension( its:ite, jts:jte ),  intent(in   )   :: &
                    iveg_ag_in                    ! combined (aggregated) dominant vegetation category   
               
               
               real(kind=kind_phys),dimension( its:ite, jts:jte ),  intent(in   )   :: &
                    firesize_avg_corr,          & ! 2d wildfire average active fire size (corrected) [m^2] read in by the host model
                    firesize_std_corr,          & ! 2d wildfire std deviation active fire size (corrected) [m^2] read in by the host model
                    frp_avg_corr,               & ! 2d wildfire average FRP (corrected) [MW] read in by the host model
                    frp_std_corr                  ! 2d wildfire std deviation FRP (corrected) [MW] read in by the host model
                    
               integer,  intent(in   )   ::     &
                    heat_flag                     ! parameter for calculate heat flux (1 - use prescribed values, 2 - use FRP)
               
               
               integer, intent(inout   ),  dimension( its:ite, jts:jte ) :: &
                    k1_zbot_out,                & ! 2d index for lower plume injection height in host model grid 
                    k2_ztop_out                   ! 2d index for upper plume injection height in host model grid  

               
               real(kind=kind_phys),dimension( its:ite, jts:jte ),  intent(out   )   :: &
                    zbot_out,                   & ! 2d lower plume injection height in host model grid        [m above surface]
                    ztop_out                      ! 2d upper plume injection height in host model grid        [m above surface]   


          !---------------------------------------------------------------
          ! Local Variables and arrays 
         

               integer ::       &
                    m1,         & ! index of the upper layer of the host model (numbered from bottom to top)
                    i,          & ! horizontal index
                    k,          & ! vertical index
                    j,          & ! j index 
                    imm           ! counter for lower and upper boundary of the heat flux. 
          
               
               real(kind=kind_phys), pointer :: &
                    ucon(:),    & ! 1D zonal wind of the host model (m/s) interpolated to theta points 
                    zcon(:),    & ! 1D cartesian height of middle of layer of the host model (meters above surface) 
                    zzcon(:),   & ! 1D cartesian height of interface of layer of the host model (meters above surface) 
                    vcon(:),    & ! 1D meridional wind of the host model (m/s) interpolated to theta points 
                    tmpcon(:),  & ! 1D Temperature of the host model (k)interpolated to theta points
                    prcon(:),   & ! 1D Pressure of the host model (Pa) interpolated to theta points
                    rvcon(:)      ! 1D water vapor mixing ratio of the host model (kg h20vapor/kg dry air)

               integer, pointer :: &
                    iveg_ag       ! combined (aggregated) dominant vegetation category 

               real(kind=kind_phys) :: &
                    burnt_area, & ! burnt area: active fire size (m2)
                    FRP           ! Fire Radiative Power (MW)

               real(kind=kind_phys), dimension(2) ::  &
                    ztopmax       ! lower and upper plume injection height at PRM grid (m above surface)

          
          call mpas_log_write(' --- enter subroutine plumerisemodel_run:')
          call mpas_log_write('     generates smoke plume injection height')

          
     !********************************************************************************
     ! Initial settings
    
     
          !--- allocate the local arrays to host model MPAS-A atmospheric conditions at 1D

               allocate(ucon(kte)  )
               allocate(zcon(kte)  )
               allocate(zzcon(kte)  )
               allocate(vcon(kte)  )
               allocate(tmpcon(kte)  )
               allocate(prcon(kte)  )
               allocate(rvcon(kte)  )
               allocate(iveg_ag) !not necessary

          !--- Define local index of the MPAS upper layer:
               
               m1=kte 

          !---- initialize several parameters as zero (only need to be done at the 1st time)
               
               call zero_plumegen_coms(ucon,zcon,zzcon,vcon,tmpcon,prcon,rvcon)

    
     !********************************************************************************
     ! Loop on host model 
            
          do j = jts,jte

               do i = its, ite

                    !---------------------------------------------------------------
                    !  if the value of active fire size is  zero, 
                    !  there is no wildfire in this grid cell,
                    !  so PRM will not be run  => cycle  
                    !  using firesize as criteria because PRM also
                    !  can run without frp. 
                    !---------------------------------------------------------------
                         
                    if (firesize_avg_corr(i, j) .lt. 1.e-6) then 
                                        
                         k1_zbot_out(i,j)=1  
                         k2_ztop_out(i,j)=1

                         cycle !go to next grid cell i

                    else

                         !call mpas_log_write('    Fire detected: run PRM ')


                         !---------------------------------------------------------------
                         ! read variables of host model  
                         !----------------------------------------------------------------
                    
                              !--- fill local 1D variables with the host model 3D variables (MPAS grid)
                              
                                   do k = kts,m1                                      ! z-direction of the host model 
                                        
                                        rvcon (k) = qv_in(i,k,j)                      ! fill with the water vapor mixing ratio (kg/kg) of the host model 
                                        zcon  (k) = hml_in(i,k,j) - z_in(i,kts,j)     ! fill with the cartesian height (meters above surface) of the termodynamic grid point  !(half-levels) above local surface--> muss full level sein, siehe aufruf htint
                                        zzcon (k) = z_in(i,k,j) - z_in(i,kts,j)       ! fill with the cartesian height (meters above surface) of the dynamic grid point  ! interface levels
                                        tmpcon(k) = temp_in(i,k,j)                    ! fill with the temperature of the host model (k) 
                                        prcon (k) = pres_in(i,k,j)                    ! fill with the pressure of the host model (Pa) . 
                                        ucon  (k) = u_in(i,k,j)                       ! fill with the zonal wind (m/s) of the host model
                                        vcon  (k) = v_in(i,k,j)                       ! fill with the meridional wind (m/s) of the host model
                                   
                                   enddo

                              !--- fill with the vegetation type combinated in the plumerismodel_driver
                              
                                   iveg_ag = iveg_ag_in(i,j)  
                                   
                         !------------------------------------------------------------------------------
                         ! get envinronmental state in the PRM grid (temp, water vapor mix ratio, ...)  
                         !----------------------------------------------------------------------------
                         
                         call get_env_condition(kte,kmt,ucon,zcon,vcon,tmpcon,prcon,rvcon)
               

                         !--------------------------------------------------------------------------------
                    
                         !Loop over the heat_flux or FRP min and max
                         do imm=1,2         
                         
                              !------------------------------------------------------------------------------
                              ! Define FireSize and FRP depending on heat_flag
                              !------------------------------------------------------------------------------
                                  
                                   !--- case heat_flag = 1 (using prescribed heat_flux)
                                   if (heat_flag == 1) then 

                                        FRP=0.0
                                        
                                        !--- case imm = 1 (lower injection height)
                                        if (imm == 1) then 

                                             burnt_area = firesize_avg_corr(i,j) - 0.674*firesize_std_corr(i,j) ! using IQR value (0.674) instead of 0.5
                                             ! --- Condition to avoid negative values ---
                                             if (burnt_area < 0.0 ) then
                                                  burnt_area = firesize_avg_corr(i,j)
                                             endif
                                        
                                        !--- case imm = 2 (upper injection height)
                                        elseif (imm == 2) then

                                             burnt_area = firesize_avg_corr(i,j) + 0.674*firesize_std_corr(i,j) ! using IQR value (0.674) instead of 0.5

                                        endif                            
                                   
                                   !--- case heat_flag = 2 (using FRP)
                                   elseif (heat_flag == 2) then

                                        !--- case imm = 1 (lower injection height)
                                        if (imm == 1) then
                                        
                                             burnt_area = firesize_avg_corr(i,j) - 0.674*firesize_std_corr(i,j) ! using IQR value (0.674) instead of 0.5
                                             FRP = frp_avg_corr(i,j) - 0.674*frp_std_corr(i,j)

                                             ! --- Condition to avoid negative values ---
                                             if (burnt_area < 0.0 .or. FRP < 0.0) then
                                                  burnt_area = firesize_avg_corr(i,j)
                                                  FRP        = frp_avg_corr(i,j)
                                             endif

                                        !--- case imm = 2 (upper injection height)
                                        elseif (imm == 2) then
                                             burnt_area = firesize_avg_corr(i,j) + 0.674*firesize_std_corr(i,j) ! using IQR value (0.674) instead of 0.5
                                             FRP = frp_avg_corr(i,j) + 0.674*frp_std_corr(i,j)
                                        endif
                                        
                                   endif


                              !------------------------------------------------------------------------------
                              ! get fire properties (burned area, plume radius, heating rates ...)  
                              !------------------------------------------------------------------------------
                                   call get_fire_properties(burnt_area, FRP, iveg_ag, heat_flag, imm)


                              !------------------------------------------------------------------------------
                              ! generate plume rise 
                              !------------------------------------------------------------------------------
                                   
                                   call makeplume (kmt,ztopmax(imm),imm)


                         enddo !end loop on imm (heat_flux min and max)        


                         !----------------------------------------------------------------------------
                         !  Print the lowest and highest plume injection height (meters above surface)
                         !----------------------------------------------------------------------------
                        
                              !DEBUG PRM - Tracking variables:
                              !call mpas_log_write('        --- lowest and highest plume height (meters above surface):')
                              !call mpas_log_write('             ztopmax(1)=$r ztopmax(2)=$r', realArgs=(/ztopmax(1),ztopmax(2)/))


                         !---------------------------------------------------------------------------
                         !  Defines the vertical layer index where the flaming phase emission will
                         !  be released in the 3d atmospheric model (MPAS grid)
                         !----------------------------------------------------------------------------
                        
                              call set_flam_vert(zzcon,ztopmax,k1,k2,zbot_plume,ztop_plume) 


                         !-------------------------------------------------------------------------------
                         ! read 1D variables from PRM run to 2D variables on PRM driver routine: 
                         !-------------------------------------------------------------------------------                                
                              
                              k1_zbot_out(i,j) = k1                   
                              k2_ztop_out(i,j) = k2                   
                              zbot_out (i,j) = zbot_plume
                              ztop_out (i,j) = ztop_plume

                    
                    endif

               enddo !end loop on i
          enddo !end loop on j

     !********************************************************************************
     ! Final settings
     !********************************************************************************

          !--- deallocate the local arrays that received MPAS atmospheric conditions at 1D

               deallocate(ucon )
               deallocate(zcon )
               deallocate(zzcon )
               deallocate(vcon )
               deallocate(tmpcon)
               deallocate(prcon)
               deallocate(rvcon)
               deallocate(iveg_ag) 



 
 call mpas_log_write(' --- end subroutine plumerisemodel_run.')
 end subroutine plumerisemodel_run
!==================================================================================================================


!=================================================================================================================
subroutine zero_plumegen_coms(ucon,zcon,zzcon, vcon,tmpcon,prcon,rvcon)

    implicit none
    real(kind=kind_phys), intent(inout) :: ucon(:),zcon(:),zzcon(:),vcon(:),tmpcon(:),prcon(:),rvcon(:)

    !call mpas_log_write('   --- enter subroutine zero_plumegen_coms:')

    w=0. ;t=0. ;td=0. ;theq=0. 
    qv=0. ;qc=0. ;qh=0. ;qi=0. ;sc=0. ;vel_p=0. ;rad_p=0. ;rad_t=0. 
    vth=0. ;vti=0. ;rho=0. ;txs=0. 
    vt3dc=0. ;vt3df=0. ;vt3dk=0. ;vt3dg=0. ;scr1=0. ;vt3dj=0. ;vt3dn=0. ;vt3do=0. 
    vt3da=0. ;scr2=0. ;vt3db=0. ;vt3dd=0. ;vt3dl=0. ;vt3dm=0. ;vt3di=0. ;vt3dh=0. ;vt3de=0. 
    wc=0. ;wt=0. ;tt=0. ;qvt=0. ;qct=0. ;qht=0. ;qit=0. ;sct=0. ;vel_t=0.
    est=0. ;qsat=0. ;qpas=0. ;qtotal=0. 
    dzm=0. ;dzt=0. ;zm=0. ;zt=0. ;vctr1=0. ;vctr2=0. 
    the=0. ;thee=0. ;rhe=0. ;sce=0. 
    ucon=0.;vcon=0.;wcon=0.;thtcon =0.;rvcon=0.;tmpcon=0.;prcon=0. 
    zcon=0.;zzcon=0.;scon=0. 
    dz=0.;dqsdz=0.;visc=0.;viscosity=0.;tstpf=0.
    advw=0.;advt=0.;advv=0.;advc=0.;advh=0.;advi=0.;cvh=0.;cvi=0.;adiabat=0.
    wbar=0.;vhrel=0.;virel=0.  
    zsurf=0.;zbase=0.;ztop=0.;area=0.;rsurf=0.;radius=0.;heating=0.
    fmoist=0.;bload=0.;dt=0.;time=0.;tdur=0.
    ztop_=0.
    n=0;nm1=0;l=0;lbase=0;mintime=0;mdur=0;maxtime=0
    dwdt_entr=0.
    kmt = 0
    vel_e(:)=0.;pe(:)=0.;te(:)=0.;qvenv(:)=0.;upe(:)=0.;vpe(:)=0.
    

    !call mpas_log_write('   --- end subroutine zero_plumegen_coms.')
   

  end subroutine zero_plumegen_coms
!=================================================================================================================


!=================================================================================================================
subroutine get_env_condition(kte,kmt,ucon,zcon,vcon,tmpcon,prcon,rvcon)

   !*******************************************************************************
    implicit none
   !*******************************************************************************

    real(kind=kind_phys), intent (in) :: & 
         ucon(:),   &                   !zonal wind (m/s) of the host model
         zcon(:),   &                   !cartesian height (meters above surface) of the termodynamic grid point  
         vcon(:),   &                   !meridional wind (m/s) of the host model 
         tmpcon(:), &                   !temperature of the host model (k) 
         prcon(:),  &                   !pressure of the host model (Pa)
         rvcon(:)                       !water vapor mixing ratio (kg/kg) of the host model 
    integer, intent (in) ::  kte        !index for the host model upper vertical grid
    integer, intent (out) :: kmt        !index for the PRM model upper vertical grid 
    
    !local:
    integer :: & 
         k,    &                        !vertical index
         nk                             !index for the host model upper vertical grid 
    
    real(kind=kind_phys) :: &
         znz                            !height of the upper termodynamic level of the host model 

    !call mpas_log_write('   --- enter subroutine get_env_condition:')

   !*******************************************************************************
   ! Define vertical grid of plume model
   !*******************************************************************************
     call set_grid()
     
     
     !this loop checks when the height of the upper vertical layer of host model is lower then PRM  
     znz=zcon(kte)  
     do k=nkp,1,-1                 !from PRM top to bottom
       if(zt(k).lt.znz) then       !compares the height of the termodynamic level of the PRM with upper termodynamic level of MPAS. Both in m above surface
          kmt=k
          exit
       endif
     enddo
    
     kmt=nkp  
     nk = kte

     !Prevents that PRM vertical grid extrapolate the host model.  
     !compares the height of the upper termodynamic level of the PRM with MPAS. Both in m above surface
       if(zt(nkp).gt.znz) then       
          stop
          call mpas_log_write('--- STOP: PRM extrapolate MPAS-A, zt(nkp)> zcon(kte):', messageType=MPAS_LOG_ERR)
       endif

   !*******************************************************************************
   ! Interpolates variables from MPAS grid to PRM grid
   !*******************************************************************************
     
     !case runing PRM accoplate to MPAS-A
     if(sound_type .ne. 1) then 
      call htint(nk,   ucon, zcon, kmt,   upe,  zt) !interpolating zonal wind
      call htint(nk,   vcon, zcon, kmt,   vpe,  zt) !interpolating meridional wind
      call htint(nk,  rvcon, zcon, kmt, qvenv,  zt) !interpolating env water vapor mixing ratio
      call htint(nk, tmpcon, zcon, kmt,    te,  zt) !interpolating temperature
      call htint(nk,  prcon, zcon, kmt,    pe,  zt) !interpolating pressure
    else 
    !case runing PRM with a sond. Why is not gonna interpolate when have the sond?
       te (:)=tmpcon(:)
       upe(:)=ucon(:)
       vpe(:)=vcon(:)
       pe (:)=prcon(:)
       qvenv(:)=rvcon(:)
    endif

   !*******************************************************************************
   ! Additional configurations: units and wind_effect
   !*******************************************************************************
     
     !Ensure that qvenv(k) does not fall below 1e-8
     do k=1,kmt
       qvenv(k)=MAX(qvenv(k),1e-8)
     enddo

     ! environmental horizontal wind (m.s-1, scalar, not the vector) computed from u and v
     do k=1,kmt 
       vel_e(k) = sqrt(upe(k)**2+vpe(k)**2)
     enddo

     !ewe - env wind effect
     if(.not. do_windeffect)  vel_e = 0.

     !converts pressure from Pa to kPa
     do k=1,kmt 
      pe(k) = pe(k)*1.e-3
     enddo

    
 !call mpas_log_write('    --- end subroutine get_env_condition.')
end subroutine get_env_condition
!=================================================================================================================

!=================================================================================================================
subroutine set_grid()

   !*******************************************************************************
    implicit none
   !*******************************************************************************
     integer :: &
         k,   & ! loop index
         mzp    ! loop PRM upper boundary

   !call mpas_log_write('   --- enter subroutine set_grid:')

   !*******************************************************************************
   ! PRM grid configuration
   !*******************************************************************************
  
    dz=100. ! set constant grid spacing of plume grid model(meters)

    mzp=nkp 
    zt(1) = zsurf 
    zm(1) = zsurf 
    zt(2) = zt(1) + 0.5*dz
    zm(2) = zm(1) + dz
    do k=3,mzp
       zt(k) = zt(k-1) + dz ! thermo and water levels
       zm(k) = zm(k-1) + dz ! dynamical levels    
    enddo
       
    do k = 1,mzp-1
       dzm(k) = 1. / (zt(k+1) - zt(k))
    enddo
    dzm(mzp)=dzm(mzp-1) 
    
    do k = 2,mzp
       dzt(k) = 1. / (zm(k) - zm(k-1))
    enddo
    dzt(1) = dzt(2) * dzt(2) / dzt(3) 
    
    !   dzm(1) = 0.5/dz
    !   dzm(2:mzp) = 1./dz

   !Print the PRM grid 
   !do k=1,mzp
     !call mpas_log_write('  zt(k)=$r', realArgs=(/zt(k)/))
   !end do


!call mpas_log_write('   --- end subroutine set_grid.')
end subroutine set_grid
!=================================================================================================================


!=================================================================================================================
subroutine htint (nzz1, vctra, eleva, nzz2, vctrb, elevb)

   !*******************************************************************************
    implicit none
   !*******************************************************************************
   
    integer,              intent(in ) :: nzz1        ! number of level from the host model (nlev)
    integer,              intent(in ) :: nzz2        ! number of level set by the global parameter nkp
    real(kind=kind_phys), intent(in ) :: vctra(nzz1) ! condition of environment from the host model 
    real(kind=kind_phys), intent(out) :: vctrb(nzz2) ! environmental conditions at PRM grid (in the plume?) -  initially set to zero by zero_plumegen_coms
    real(kind=kind_phys), intent(in ) :: eleva(nzz1) ! vertical grid of the host model
    real(kind=kind_phys), intent(in ) :: elevb(nzz2) ! vertical grid set by set_grid()
    
    integer  :: l  ! counter  
    integer  :: k  ! counter
    integer  :: kk ! counter for print(log.write) 
    real(kind=kind_phys) :: wt !weight for grid interpolation

    !call mpas_log_write('   --- enter subroutine htint: grid interpolation')

   !*******************************************************************************
   !Grid interpolation
   !*******************************************************************************
   
     l=1

     do k=1,nzz2
       do
          if ( (elevb(k) <  eleva(1)) .or. &
               ((elevb(k) >= eleva(l)) .and. (elevb(k) <= eleva(l+1))) ) then
             wt       = (elevb(k)-eleva(l))/(eleva(l+1)-eleva(l))
             vctrb(k) = vctra(l)+(vctra(l+1)-vctra(l))*wt
             exit
          else if ( elevb(k) >  eleva(nzz1))  then
             wt       = (elevb(k)-eleva(nzz1))/(eleva(nzz1-1)-eleva(nzz1))
             vctrb(k) = vctra(nzz1)+(vctra(nzz1-1)-vctra(nzz1))*wt
             exit
          end if
          
          l=l+1
          if(l == nzz1) then
             ! PRINT *,'htint:nzz1',nzz1 !commented out for vectorisation 
             ! DO kk=1,l
             !   PRINT*,'kk,eleva(kk),elevb(kk)',eleva(kk),elevb(kk)
             ! END DO
             ! STOP 'htint'
             EXIT ! NP: Exit instead of Stop
          end if
       end do
     end do

!call mpas_log_write('   --- end subroutine htint: grid interpolation.')
end subroutine htint
!=================================================================================================================


!=================================================================================================================
subroutine get_fire_properties(burnt_area, FRP, iveg_ag, heat_flag, imm)

   !*******************************************************************************
    implicit none
   !*******************************************************************************
   
    real(kind=kind_phys), intent(in) :: &
         burnt_area,     &  ! burnt area read in by the host model
         FRP                ! Fire Radiative Power read in by the host model
       
    integer, intent(in) :: iveg_ag, heat_flag, imm  
    

    integer :: &
         moist,   & ! fuel moisture
         i,       & ! loop counter ntime
         icount     ! loop counter number of seconds in the burn
    
    real(kind=kind_phys) :: &
         bfract,        & ! combustion factor
         effload,       & ! patchy burning NP: not used at all
         heat,          & ! joules/kg - tropical forest (mt) NP: there are other alternatives in the comments but only this one is used
         hinc             ! heat increase for a 5 minute ramp up in 4 steps
    
    real(kind=kind_phys), parameter :: beta = 0.88  !ref.: Paugam  et al., 2015 
     
    real(kind=kind_phys),  dimension(2,4) :: heat_flux
         data heat_flux/  &
         !---------------------------------------------------------------------
         !  heat flux         !IGBP Land Cover      ! 
         !   kW/m^2           !Legend and           ! reference
         ! min   max          !description          ! 
         !--------------------------------------------------------------------
           30.0,  80.0,   &! iveg_ag=1:Tropical Forest         ! igbp 2
           30.0,  80.0,   &! iveg_ag=2:Boreal forest           ! igbp 1 & 3,4,5
            4.4,  23.0,   &! iveg_ag=3:cerrado/woody savanna   | igbp  6 thru 9, and 11
            3.3,   3.3    /! iveg_ag=4:Grassland/cropland      ! igbp 10, 12, 14
         !--------------------------------------------------------------------

     !call mpas_log_write('   --- enter subroutine get_fire_properties:')

   !*******************************************************************************
   ! Fire data inputs
   !*******************************************************************************

    !-- fire area at the surface
    area = burnt_area   ! area of burn, m^2
    

    !-- heat flux for the vegetation type
    if ( heat_flag == 1) then

         heat_fluxW = heat_flux(imm,iveg_ag) * 1000. !take the correspondent value in the table and convert to W/m^2
      
     elseif ( heat_flag == 2) then
         
         ! "beta" factor converts FRP to convective energy:
         !Check with Saulo
         heat_fluxW = beta * (FRP * (10.0**6) / area) / 0.55 ! in W/m^2 !Jaqueline Pereira (add the conversion *10^6 to convert FRP from MW to W)
    endif

    

    !-- specific heat joules per kg of fuel consumed  
    !Why use the same value for the all vegetation since we have differents values?
    !heat = 21.e6     !-                  
    !heat = 15.5e6    !- joules/kg - cerrado/savannah
    heat = 19.3e6     !- joules/kg - tropical forest (mt)


    !-- others fires parameters
    mdur = 53        ! duration of burn, minutes
    !bload = 10.      ! total loading, kg/m**2 : not used at all 
    moist = 10       ! fuel moisture, %. average fuel moisture,percent dry
    maxtime =mdur-1  ! model time, min


   !*******************************************************************************
   ! Fix up inputs
   !*******************************************************************************

    if (MOD (maxtime, 2) .ne.0) maxtime = maxtime+1  !make maxtime even
    maxtime = maxtime * 60    ! and put in seconds

    rsurf = sqrt (area / pi)  !- entrainment surface radius (m)

    fmoist   = moist / 100.   !- fuel moisture fraction (%)


    !------------------------------------------------------------------------------
    ! calculate the energy flux and water content at lboundary.
    !----------------------------------------------------------------------------

     !fills heating() on a minute basis. 
     !could ask for a file at this po in the program. 
     !whatever is input has to be adjusted to a one minute timescale.

     do i = 1, ntime            !- make sure of energy release
       heating (i) = 0.0001     !- avoid possible divide by 0
     enddo

     tdur = mdur * 60.          !- number of seconds in the burn

     !bfract = 1.                !- combustion factor: not used at all

     !effload = bload * bfract   !- patchy burning: not used at all


    !------------------------------------------------------------------------------
    ! spread the burning evenly over the interval of the burn duration
    !----------------------------------------------------------------------------
     ! for the first few minutes (5min) the convective energy is lower for stability

     icount = 1  
     if(mdur > ntime) STOP 'Increase time duration (ntime) in min - see file "mpas_chemistry_plumerisemodel_run.F90"'

     do while (icount.le.mdur)                             
          heating (icount) = heat_fluxW  * 0.55     ! W/m**2 (0.55 converts to convective energy)
          icount = icount + 1  
     enddo
          
    !Ramp for 5 minutes
    !The convective energy increases linearly from 0.1 to heat_fluxW*0.55 in the first 5min. 
    hinc = heating (1) / 4.  
    heating (1) = 0.1              !convective energy = 0.1
    heating (2) = hinc             !convective energy = (heat_fluxW  * 0.55)*0.25
    heating (3) = 2. * hinc        !convective energy = (heat_fluxW  * 0.55)*0.50
    heating (4) = 3. * hinc        !convective energy = (heat_fluxW  * 0.55)*0.75


    !call mpas_log_write('   --- end subroutine get_fire_properties.')

 return 
end subroutine get_fire_properties
!=================================================================================================================


!=================================================================================================================
subroutine makeplume(kmt,ztopmax,imm)

     ! *********************************************************************
     !
     !    EQUATION SOURCE--Kessler Met.Monograph No. 32 V.10 (K)
     !    Alan Weinstein, JAS V.27 pp 246-255. (W),
     !    Ogura and Takahashi, Monthly Weather Review V.99,pp895-911 (OT)
     !    Roger Pielke,Mesoscale Meteorological Modeling,Academic Press,1984
     !    Originally developed by: Don Latham (USFS)
     !
     !
     ! ------------------------ VARIABLE ID -----------------------------
     !
     !     DT=COMPUTING TIME INCREMENT (SEC)
     !     DZ=VERTICAL INCREMENT (M)
     !     LBASE=LEVEL ,CLOUD BASE
     !
     !     CONSTANTS:
     !       G = GRAVITATIONAL ACCELERATION 9.80796 (M/SEC/SEC).
     !       R = DRY AIR GAS CONSTANT (287.04E6 JOULE/KG/DEG K)
     !       CP = SPECifIC HT. (1004 JOULE/KG/DEG K)
     !       HEATCOND = HEAT OF CONDENSATION (2.5E6 JOULE/KG)
     !       HEATFUS = HEAT OF FUSION (3.336E5 JOULE/KG)
     !       HEATSUBL = HEAT OF SUBLIMATION (2.83396E6 JOULE/KG)
     !       EPS = RATIO OF MOL.WT. OF WATER VAPOR TO THAT OF DRY AIR (0.622)
     !       DES = DifFERENCE BETWEEN VAPOR PRESSURE OVER WATER AND ICE (MB)
     !       TFREEZE = FREEZING TEMPERATURE (K)
     !
     !  
     !     PARCEL VALUES:
     !       T = TEMPERATURE (K)
     !       TXS = TEMPERATURE EXCESS (K)
     !       QH = HYDROMETEOR WATER CONTENT (G/G DRY AIR)
     !       QHI = HYDROMETEOR ICE CONTENT (G/G DRY AIR)
     !       QC = WATER CONTENT (G/G DRY AIR)
     !       QVAP = WATER VAPOR MIXING RATIO (G/G DRY AIR)
     !       QSAT = SATURATION MIXING RATIO (G/G DRY AIR)
     !       RHO = DRY AIR DENSITY (G/M**3) MASSES = RHO*Q'S IN G/M**3
     !       ES = SATURATION VAPOR PRESSURE (kPa)
     !
     !     ENVIRONMENT VALUES:
     !       TE = TEMPERATURE (K)
     !       PE = PRESSURE (kPa)
     !       QVENV = WATER VAPOR (G/G)
     !       RHE = RELATIVE HUMIDITY FRACTION (e/esat)
     !       DNE = dry air density (kg/m^3)
     ! 
     !     HEAT VALUES:
     !       HEATING = HEAT OUTPUT OF FIRE (WATTS/M**2)
     !       MDUR = DURATION OF BURN, MINUTES
     !
     !       W = VERTICAL VELOCITY (M/S)
     !       RADIUS=ENTRAINMENT RADIUS (FCN OF Z)
     !       RSURF = ENTRAINMENT RADIUS AT GROUND (SIMPLE PLUME, TURNER)
     !       ALPHA = ENTRAINMENT CONSTANT
     !       MAXTIME = TERMINATION TIME (MIN)
     !
     !
     !**********************************************************************



     !*******************************************************************************
      implicit none
     !*******************************************************************************

      integer, intent (in) :: &
          kmt,        & ! set to nkp = 200 in get_env_condition
          imm           !
      real(kind=kind_phys), intent(out) :: ztopmax
      
      character (len=10) :: varn

      integer :: &
          iconv,      & !
          itime,      & !
          k,          & !
          kk,         & !
          kkmax,      & !
          deltak,     & !
          nrectotal,  & ! 
          i_micro,    & !
          n_sub_step, & !
          ilastprint
     
      real(kind=kind_phys) :: &
          wmax,     & !
          rmaxtime, & !
          es,       & !
          heat,     & !
          dt_save     !
      real(kind=kind_phys), parameter :: vc = 5.107387  
      real(kind=kind_phys), parameter :: &
          g = 9.80796, & ! LM: from global parameters
          r = 287.04,  & ! LM: from global parameters
          cp = 1004.,  & ! LM: from global parameters
          eps = 0.622, & !
          tmelt = 273.3  ! LM: from global parameters
      real(kind=kind_phys), parameter :: &
          heatsubl = 2.834e6, & !
          heatcond = 2.501e6    !
      
      ! -----------------------------------------------------------
      ! some constants 
      !-------------------------------------------------------------
      !vc = 38.3/(xno**.125) mean volume fallspeed eqn. (K)
      !xno=10.0e06 median volume diameter raindrop (K table 4)

      
      !call mpas_log_write('   --- enter subroutine makeplume:')     

     !*******************************************************************************
     ! Problem setup and initial conditions
     !*******************************************************************************
      tstpf = 2.0     !- timestep factor
      nrectotal=kmt   !-kmt (niko) or 150 (saulo)? 

      mintime = 1  
      ztopmax = 0. 
      ztop    = 0. 
      time = 0.  
      dt = 1.
      wmax = 1. 
      kkmax   = 10
      deltaK  = 20
      ilastprint=0
      l       = 1   ! l initialization

      viscosity = 500.!- viscosity constant (original value: 0.001)


      ! -----------------------------------------------------------
      ! initialization
      !-------------------------------------------------------------
           
          call initial(kmt)
             
     
        !*******************************************************************************
        ! Model Evolution 
        !*******************************************************************************
        
          
            ! -----------------------------------------------------------
            ! Begins  Time Loop 
            !-------------------------------------------------------------
            
            
            rmaxtime = float(maxtime)

            !Beginning of time loop
            DO while (time.le.rmaxtime)  !while time(starts as zero) lower or equal 53min -1 (duration of burn)
                
                !-- set model top integration
                !Defines the top of the vertical layer to integration as the minimum value between
                !kmt=200 and kkmax+deltak that starts as 10+20=30, which means that will initiate
                !the integration from surface to the level k=30. This means that the model will 
                !integrate just the cloud instead the whole vertical grid. 
                nm1 = MIN(kmt, kkmax + deltak)

                
                !-- set timestep in seconds
                !calculates the timestep based on spatial parameters and a timestep factor.
                !The calculated dt is capped at a maximum of 5sec to prevent excessively large timesteps. 
                dt = MIN(5.,((zm(2)-zm(1)) / (tstpf * wmax))) 
                
                !why niko removed the !-PNNL recomm section here, which is in brams(except the new dt)?
                !-PNNL recomm
                        !dt = min(dt,0.8 * 3.14159 * radius(2) / (2.0 * abs(vel_p(2) - vel_e(2))))
                !-PNNL recomm
                        !IF(dt < 0.25 .or. maxval(T(:))<0.) THEN
                        !     ztop = 0.
                        !     ztop_(mintime) = ztop
                        !     ztopmax = 0.
                        !     kkmax	= 1 
                        !     print*,"PRM_UNSTABLE=",dt
                        !     !call flush(6)
                        !     EXIT
                        !ENDIF

                !-- elapsed time, seconds
                !updates the time
                time = time+dt
                
                !-- elapsed time, minutes
                !Converts the elapsed time from seconds to minutes. 
                !The integer division by 60 discards any fractions, and adding 1 ensures mintime starts at 1.                                      
                mintime = 1 + INT (time) / 60     
                wmax = 1.  !no zeroes allowed.
                

               ! -----------------------------------------------------------
               ! Begins Space Loop 
               !-------------------------------------------------------------

                !-- zerout all model tendencies
                call tend0_plumerise ! sets tendency variables to zero
            

                !-- surface boundary conditions (k=1)
                l=1
                !wbar=w(1) !just in brams. i need this here?
                call lbound_mtt() ! sets lower bound entries. 


                !-- dynamics for the level k>1 
                !-- W advection 
                call vel_advectc_plumerise(nm1,wc,wt,rho,dzm)


                !-- scalars advection (advection for scalars: qv, qi, temperature, etc)
                !case original scheme
                if(.not. use_bott_adv) then
                        call scl_advectc_plumerise (nm1)
                else
                !case Bott's advection scheme
                        call scl_advectc_plumerise2(nm1)
                endif


                !-- scalars entrainment, adiabatic
                call scl_misc(nm1) 


                !-- scalars dinamic entrainment
                call  scl_dyn_entrain(nm1)


                !-- gravity wave damping using Rayleigh friction layer fot T
                call damp_grav_wave(1,nm1,deltak,zt,w,t,tt,qv,qh,qi,qc,te,pe,qvenv, &
                                    vel_p,vel_t,vel_e) 
                
                
                !-- microphysics
                if(do_microphysics) then     
                                  
                    dt_save=dt
                    n_sub_step=3
                    dt=dt/float(n_sub_step)
                    
                    do i_micro=1,n_sub_step
                         !-- sedim
                         call fallpart(nm1)
                         do l=2,nm1 !the limits here are a little different from the original: l=2,nm1-1.  ?
                              es      = 0.1*esat(t(l)) !blob saturation vapor pressure, em kPa
                              qsat(l) = (eps * es) / (pe(l) - es)  !blob saturation lwc g/g dry air
                              est (l) = es  
                              rho (l) = 3483.8 * pe (l) / t (l) ! AIR PARCEL DENSITY , G/M**3
                         enddo
                         
                         do l=2,nm1-1
                              wbar    = 0.5*(w(l)+w(l-1))
                              dqsdz(l) = (qsat(l+1) - qsat(l-1)) / (zt(l+1) -zt(l-1))
                              call waterbal  
                         enddo
                    enddo

                    dt=dt_save
                endif 


                !-- W-viscosity for stability 
                call visc_W(nm1,deltak,kmt)


                !-- update scalars
                call update_plumerise(nm1,'S')
                call hadvance_plumerise(1,nm1,wc,wt,w,mintime)


                !-- Buoyancy
                call buoyancy_plumerise(nm1, t, te, qv, qvenv, qh, qi, qc, wt, scr1)


                !-- Entrainment 
                call entrainment(nm1,w,wt,radius,alpha)


                !-- update W
                call update_plumerise(nm1,'W')
                call hadvance_plumerise(2,nm1,wc,wt,w,mintime) 


                !-- misc . i dont understand this part yet. ?
                do k=2,nm1
                    !    pe esta em kpa  - esat do rams esta em mbar = 100 Pa = 0.1 kpa
                    !    I need to change according to mpas?
                    es       = 0.1*esat(t(k)) !blob saturation vapor pressure, em kPa
                    !    rotina do plumegen calcula em kPa
                    qsat(k) = (eps * es) / (pe(k) - es)  !blob saturation lwc g/g dry air
                    est (k) = es  
                    txs (k) = t(k) - te(k)
                    rho (k) = 3483.8 * pe (k) / t (k) ! air parcel density , g/m**3
                                                      ! no pressure diff with radius

                    if((ABS(wc(k))).gt.wmax) wmax = ABS(wc(k)) ! keep wmax largest w
                enddo

                ! Gravity wave damping using Rayleigh friction layer for W
                call damp_grav_wave(2,nm1,deltak,zt,w,t,tt,qv,qh,qi,qc,te,pe,qvenv, &
                                   vel_p,vel_t,vel_e) !again it's not passing dt and zm.?

                !- update radius
                do k=2,nm1
                    radius(k) = rad_p(k)
                enddo

                !-- try to find the plume top (above surface height)
          
                kk = 1
                do while (w (kk) .gt. 1.)  
                    kk = kk + 1  
                    ztop =  zm(kk)
                enddo
               
                ztop_(mintime) = ztop
                ztopmax = MAX (ztop, ztopmax) 
                !updates the top for integration 
                kkmax   = MAX (kk  , kkmax  )
                
                !-- case the solution is going to a stationary phase: 
                if(mintime > 40) then
                     if( abs(ztop_(mintime)-ztop_(mintime-10)) < dz ) exit
                endif
             

            ENDDO !go to next timestep


!call mpas_log_write('    --- end subroutine makeplume.') 
end subroutine makeplume
!=================================================================================================================


!=================================================================================================================
subroutine initial(kmt)

     !*******************************************************************************
      implicit none
     !*******************************************************************************
      integer, intent(in) :: kmt ! set to nkp = 200 in get_env_condition
      integer :: k  ! loop counter from 1 to N=kmt
      real(kind=kind_phys) :: es     ! blob saturation vapor pressure, em kPa
      real(kind=kind_phys) :: x,y


      !call mpas_log_write('   --- enter subroutine initial:')

     !*******************************************************************************
     ! initialize temperature structure,to the end of equal spaced sounding
     !*******************************************************************************
      n=kmt
         do k = 1, n             
           txs (k) = 0.0  
           w (k) = 0.0             
           t (k) = te(k) !blob set to environment
           wc(k) = 0.0
           wt(k) = 0.0
           qv(k) = qvenv (k)                      
           vth(k) = 0.      !initial rain velocity = 0                       
           vti(k) = 0.      !initial ice  velocity = 0                       
           qh(k) = 0.       !no rain                 
           qi(k) = 0.       !no ice                  
           qc(k) = 0.       !no cloud drops                      
           es    = 0.1*esat(t(k)) !blob saturation vapor pressure, em kPa
           est  (k) = es  
           qsat (k) = (.622 * es) / (pe (k) - es) !saturation lwc g/g
           rho  (k) = 3483.8 * pe (k) / t (k)     !dry air density g/m**3     
           vel_p(k) = 0.
           rad_p(k) = 0.
         enddo

     !*******************************************************************************
     ! initialize the entrainment radius
     !*******************************************************************************

      ! alpha
      alpha(:) = alpha_mtt

      ! -----------------------------------------------------------
      ! Tarhish and Romps plume option 
      !-------------------------------------------------------------
         if (.not. use_zv_mtt) then 
            do k=1,n
               x = zm(k)/rsurf
             if(zt(k) < 2.*rsurf) then 
               y = 0.02369016 + (0.5199647 - 0.02369016)/(1 + (x/0.2040715)**0.8115766)
               alpha(k) = 1.5*max(y,alpha_mtt)
             endif
             !print*,'alpha=',real(zm(k),4),real(alpha(k),4) ,real(x,4),real(y,4)
            enddo
         endif
    
      ! -----------------------------------------------------------
      ! Turner plume option 
      !-------------------------------------------------------------

          radius(1) = rsurf
          rad_p (1) = rsurf 
          do k=2,n
               radius(k) = radius(k-1)+(6./5.)*alpha(k)*(zt(k)-zt(k-1)) !aqui pode deixar o alpha dependente k mesmo sendo constante? 
               rad_p(k)  = radius(k)
          enddo

     !*******************************************************************************
     ! initialize the viscosity
     !*******************************************************************************
      !values of visc go from 500 to ? with k=1 to k=200
      visc (1) = viscosity
      do k=2,n
          visc (k) = max(1.e-3,visc(k-1) - 1.* viscosity/real(nkp, kind_phys))
      enddo


     !*******************************************************************************
     ! initialize lower boundary conditions
     !*******************************************************************************
      call lbound_mtt()

     !Initialize non-temperature variables? rotina nao presente no modelo niko.
     !Check with Saulo
     

 !call mpas_log_write('    --- end subroutine initial.')
end subroutine initial
!=================================================================================================================


!=================================================================================================================
subroutine lbound_mtt()

     ! ************************************************************************** 
     !          BOUNDARY CONDITIONS AT ZSURF FOR PLUME AND CLOUD 
     ! source of equations: J.S. Turner Buoyancy Effects in Fluids
     !                      Cambridge U.P. 1973 p.172,
     !                      G.A. Briggs Plume Rise, USAtomic Energy Commissio
     !                      TID-25075, 1969, P.28
     ! fundamentally a point source below ground. at surface, this produces
     ! a velocity w(1) and temperature T(1) which vary with time. There is
     ! also a water load which will first saturate, then remainder go into
     ! QC(1).
     ! EFLUX = energy flux at ground,watt/m**2 for the last DT
     ! ************************************************************************** 
     

     !*******************************************************************************
      implicit none
     !*******************************************************************************
      real(kind=kind_phys), parameter :: &
         g = 9.80796, & ! from global parameters
         r = 287.04,  & ! from global parameters
         cp = 1004.6, & ! from global parameters
         eps = 0.622    ! 

      real(kind=kind_phys), parameter :: &
         pi = 3.14159,    & ! from global parameters
         e1 = 1./3.,      & !
         e2 = 5./3.         !

      real(kind=kind_phys) :: &
         es,      & ! blob saturation vapor pressure, em kPa
         eflux,   & ! energy flux at lboundary computed in burn(.,.)
         water,   & ! water set free by fire at lboundary and then computed in burn(.,.)
         pres,    & ! pressure in N/m**2
         c1,      & ! 5. / (6. * entrainment constant)
         c2,      & ! 0.9 * entrainment constant
         f,       & ! buoyancy flux
         zv,      & ! virtual boundary height
         denscor, & ! density correction
         xwater     ! firewater mixing ratio

         !esat_pr? nao tava no modelo do niko. Check with Saulo

      !call mpas_log_write('    --- enter subroutine lbound_mtt:')

     !*******************************************************************************
     ! calculate parameters at boundary from a virtual buoyancy point source
     !*******************************************************************************

      ! -----------------------------------------------------------
      !calculates the energy flux and water content at lboundary
      !-------------------------------------------------------------

          qh (1) = qh (2)   !soak up hydrometeors
          qi (1) = qi (2)              
          qc (1) = 0.       !no cloud here

          call burn (eflux, water) 
       
      ! -----------------------------------------------------------
      !calculates boundary conditions
      !-------------------------------------------------------------
       
          !-- need pressure in N/m**2
               pres = pe (1) * 1000.
     
          !-- parameters for the lower boundary equations
               c1 =  5. / (6. * alpha(1))  
               c2 = (9. /10.) * alpha(1)  

          !-- virtual boundary height 
               
               !(Latham, 1994)
               zv = -c1 * rsurf
               
               !-- Tarshish & Romps (2022, 10.1175/JAS-D-21-0096.1)
               if(.not. use_zv_mtt) zv = -rsurf

          !-- buoyancy flux  (@ z = 0)
               f = eflux / (pres * cpd)  
               f = grav * rd * f * (area/pi)  

          !--  boundary condition for vertical velocity (@ z = 0)
               w (1) = c1 * ( (c2 * f) **e1) / (-zv)**e1  
               

          !-- density correction  (@ z = 0)
               denscor = c1 * f / grav / (c2 * f) **e1 / (-zv)**e2   

          !-- temperature of virtual plume (@ z = 0)
               t (1) = te (1) / (1. - denscor)   

          !-- boundary condition for vertical velocity, plume radius, etc  (@ z = 0)
               
               wc   (1) = w(1)
               vel_p(1) = 0.
               rad_p(1) = rsurf
               sc   (1) = 1000. !--gas/particle (g/g). This is correct? !pode ignorar

               vth (1) = - 4.  
               vti (1) = - 3.  
               txs (1) = t (1) - te (1)  
               visc (1) = viscosity  


          !-- air density at level 1, g/m**3
               rho (1) = 3483.8 * pe (1) / t (1)  
    
          !-- firewater mixing ratio
               xwater = water / (w (1) * dt * rho (1) ) 

          !-- boundary condition for water vapor 
               
               qv (1) = xwater + qvenv (1)  !plus what's already there 
               
               
               !--  pe esta em kPa  - esat esta em mbar = 100 Pa = 0.1 kPa
               es       = 0.1*esat(t(1)) !blob saturation vapor pressure, em kPa
               est  (1) = es                                  
               qsat (1) = (eps * es) / (pe (1) - es)   !blob saturation lwc g/g dry air

               if (qv (1) .gt. qsat (1) ) then  
                    qc (1) = qv   (1) - qsat (1) + qc (1)  !remainder goes into cloud drops
                    qv (1) = qsat (1)  
               endif

          !-----------------------------------------------------------------
          !microphyscis species
          !-----------------------------------------------------------------
           call waterbal

  !call mpas_log_write('    --- end subroutine lbound_mtt:')
 return
end subroutine lbound_mtt
!=================================================================================================================


!=================================================================================================================
subroutine burn(eflux, water)

     !*******************************************************************************
     ! The emission factor for water is 0.5. The water produced, in kg,
     ! is then  fuel mass*0.5 + (moist/100)*mass per square meter.
     ! The fire burns for DT out of TDUR seconds, the total amount of
     ! fuel burned is AREA*BLOAD*(DT/TDUR) kg. this amount of fuel is
     ! considered to be spread over area AREA and so the mass burned per
     ! unit area is BLOAD*(DT/TDUR), and the rate is BLOAD/TDUR.
     !******************************************************************************
     

     !*******************************************************************************
      implicit none
     !*******************************************************************************
      real(kind=kind_phys), parameter  :: &
          !heat = 21.E6  , & !Joules/kg
          !heat = 15.5E6 , & !Joules/kg - cerrado
          heat = 19.3e6      !Joules/kg - floresta em Alta Floresta (MT)
      real(kind=kind_phys), intent(out) :: &
          eflux, & ! energy flux at lboundary
          water    ! water set free by fire at lboundary

      !call mpas_log_write('    --- enter subroutine burn:')

     !*******************************************************************************
     ! calculates the energy flux and water content at lboundary
     !*******************************************************************************
     
          if (time.gt.tdur) then !is the burn over?   
               eflux = 0.000001    !prevent a potential divide by zero
               water = 0.  
               return  
          else                   
               eflux = heating (mintime)                              ! Watts/m**2                                                   
               !WATER = EFLUX * (DT / HEAT) * (0.5 + FMOIST)          ! kg/m**2 !when use this equation?
               water = eflux * (dt / heat) * (0.5 + fmoist) /0.55     ! kg/m**2 
               water = water * 1000.                                  ! g/m**2
          endif
 
 !call mpas_log_write('    --- end subroutine burn:')

 return 
end subroutine burn
!=================================================================================================================

!=================================================================================================================
subroutine waterbal

 implicit none 

     if (qc (l) .le. 1.0e-10) qc (l) = 0.  !defeat underflow problem 
     if (qh (l) .le. 1.0e-10) qh (l) = 0.  
     if (qi (l) .le. 1.0e-10) qi (l) = 0.  

     !call mpas_log_write('    --- enter subroutine waterbal:')

     call evaporate    !vapor to cloud,cloud to vapor  

     call sublimate    !vapor to ice  

     call glaciate     !rain to ice 

     call melt         !ice to rain

     !if(ak1 > 0. .or. ak2 > 0.) &                      !I dont understand?
     call convert () !(auto)conversion and accretion 
     !call convert2 () !(auto)conversion and accretion  !I dont understand? I need to bring this subroutine here?

    !call mpas_log_write('    --- end subroutine waterbal.')

 return 
end subroutine waterbal
!=================================================================================================================

!=================================================================================================================
subroutine evaporate
     
     !- evaporates cloud,rain and ice to saturation

     !*****************************************************************
      implicit none
     !*****************************************************************
      real(kind=kind_phys), parameter :: &
         herc = 5.44e-4,  & ! evaporation constant, HERC = 1.93*1.E-6*XN035, XNO=10.0E06
         cp = 1.004,      & ! from global parameters
         heatcond = 2.5e3   ! from global parameters
      real(kind=kind_phys), parameter :: &
         heatsubl = 2834.   !
      real(kind=kind_phys), parameter :: &
         frc = heatcond / cp, & ! 
         src = heatsubl / cp    !
      real(kind=kind_phys) :: &
         evhdt,    & !
         evidt,    & !
         evrate,   & !
         evap,     & !
         sd,       & !
         quant,    & !
         dividend, & !
         divisor,  & !
         devidt      !

      !call mpas_log_write('    --- enter subroutine evaporate:')

     !*****************************************************************
     ! calculations for vapor deficit and eavaporation rate
     !*****************************************************************
         
          !vapor deficit
          sd = qsat (l) - qv (l) 

          !if air is saturated, no evaporation or condesation occur 
          if (sd.eq.0.0)  return  
          !if (abs(SD).lt.1.e-7)  RETURN  

          evhdt = 0.  
          evidt = 0. 
          
          !evrate =0.; evap=0.; sd=0.0; quant=0.0; dividend=0.0; divisor=0.0; devidt=0.0   !why this is commented?                           
          
          ! evaporation rate (Kessler 8.32)
          evrate = ABS (wbar * dqsdz(l))   

          ! total amount of water available to evaporate in a given dt
          evap = evrate * dt   

     !*****************************************************************
     ! Air supersaturate: condensation process starts
     !*****************************************************************
          
         ! case sd negative, qv>qsat, so condense.
         IF (sd.le.0.0) then   
       
               !case evap>=|sd|, all deficit is used to condense:
               if (evap.ge. ABS (sd) ) then                                   
                    qc (l) = qc  (l) - sd         !deficit,remember? means actually qc + |sd|
                    qv (l) = qsat(l)              !set the vapor to saturation  
                    t  (l) = t   (l) - sd * frc   !heat gained through condensation per gram of dry air. deficit,remember? means actually t + |sd|*frc
                    return  
                else
                !case evap<|sd|, just the evap is used to condense:   
                    qc (l) = qc (l) + evap         !get what we can in dt 
                    qv (l) = qv (l) - evap         !remove it from the vapor
                    t  (l) = t  (l) + evap * frc   !get some heat
                    return
               endif 
               
         ELSE
                  
           !*****************************************************************
           ! Air unsaturate: evaporation process starts to reach saturation
           !*****************************************************************

            ! case sd is positive, qv<qsat, so evaporate.
            ! saturate if possible. 
            ! use everything in order: cloud, rain, ice.
            ! in case of cloud evaporation is not enough to saturate the air, move to rain and then ice. 

            !-------------------------------------------------------
            !evaporation of cloud water
            !-------------------------------------------------------
           
               !case evap <= qc, so part of cloud droplets will evaporate
               if (evap.le.qc (l) ) then           !enough cloud to last DT  
                
                !evap is enough to saturate the air:
                if (sd.le.evap) then           !enough time to saturate
                    qc (l) = qc (l) - sd       !remove cloud                                          
                    qv (l) = qsat (l)          !saturate
                    t (l) = t (l) - sd * frc   !cool the parcel                                          
                    return   !done                    
                
                !evap not enough to saturate the air:
                else   !not enough time                                 
                    sd = sd-evap               !use what there is. update the deficit. the deficit will become lower because of the evaporation
                    qv (l) = qv (l) + evap     !add vapor
                    t (l) = t (l) - evap * frc !lose heat, cool the parcel.
                    qc (l) = qc (l) - evap     !lose cloud
                    !go on to rain.                                      
                endif
           
               !case evap >= qc, all clouds droplets will evaporate. not enough cloud to last dt 
               else 

                !case sd<=qc, all the cloud droplets will evaporate and saturate the air
                if (sd.le.qc (l) ) then                                     
                    qv (l) = qsat (l)              !use it
                    qc (l) = qc (l) - sd  
                    t  (l) = t (l) - sd * frc  
                    return                             
                !case sd>qc, sll the cloud droplets will evaporate but not saturate the air
                else            
                    sd = sd-qc (l)  
                    qv (l) = qv (l) + qc (l)  
                    t  (l) = t (l) - qc (l) * frc         
                    qc (l) = 0.0                         !all gone                                  
                endif       !on to rain                          
               endif          !finished with cloud


           !-------------------------------------------------------
           !evaporation of rain
           !-------------------------------------------------------
           
            !  if all cloud droplets evaporate but still not saturate the air, 
            !  so try to use some rain.
            !  this is tricky, because we only have time DT to evaporate. 
            !  if there is enough rain, we can evaporate it for dt. 
            !  ice can also sublimate  at the same time. 
            !  there is a compromise here.....use rain first, then  ice. 
            !  saturation may not be possible in one DT time.
            !  rain evaporation rate (W12),(OT25),(K Table 4). evaporate rain first
            !  sd is still positive or we wouldn't be here.

               !check if there is rain
               if (qh (l) .gt.1.e-10) then 
               
                    !srf-25082005 
               
                    ! amount of water to saturate the air after evaporation of cloud droplets
                    !  QUANT = (QC (L) + QV (L) - QSAT (L) ) * RHO (L)   !g/m**3 pq o sinal aqui ta invertido?
                    quant = ( qsat (l)-qc (l) - qv (l)  ) * rho (l)   !g/m**3
               
                    !amount of rain available to evaporate in time DT :
                    evhdt = (dt * herc * (quant) * (qh (l) * rho (l) ) **.65) / rho (l)                
               

                    if (evhdt.le.qh (l) ) then           !enough rain to last dt
                     
                     !case the rain available to evaporate is enough to saturate the air:
                     if (sd.le.evhdt) then           !enough time to saturate      
                         qh (l) = qh (l) - sd       !remove rain      
                         qv (l) = qsat (l)          !saturate     
                         t (l) = t (l) - sd * frc   !cool the parcel          
                         return              !done                      
                     else                               !not enough time
                         sd = sd-evhdt                 !use what there is
                         qv (l) = qv (l) + evhdt       !add vapor
                         t (l) = t (l) - evhdt * frc   !lose heat
                         qh (l) = qh (l) - evhdt       !lose rain   
                     endif                   !go on to ice.
               
                    !case  evhdt>=qh, all rain will evaporate
                    else 
                               
                    !case all rain will evaporate and saturate the air:
                     if (sd.le.qh (l) ) then           !but there is enough to saturate
                         qv (l) = qsat (l)             !use it
                         qh (l) = qh (l) - sd  
                         t (l) = t (l) - sd * frc  
                         return                            
                
                     !all rain will evaporate but not saturate the air:
                     else                              !not enough to sat
                         sd = sd-qh (l)  
                         qv (l) = qv (l) + qh (l)  
                         t (l) = t (l) - qh (l) * frc    
                         qh (l) = 0.0                   !all gone                                          
                     endif                              !on to ice                                       
                    
                    endif                                !finished with rain
             endif 


           !-------------------------------------------------------
           !evaporation of ice
           !-------------------------------------------------------
           
            ! now for ice
            ! equation from (OT); correction factors for units applied

            if (qi (l) .le.1.e-10) return            !no ice there

            !amount of ice available to evaporate (rate)
            dividend = ( (1.e6 / rho (l) ) **0.475) * (sd / qsat (l) - 1) * (qi (l) **0.525) * 1.13
            divisor = 7.e5 + 4.1e6 / (10. * est (l) )  
            devidt = - cvi(l) * dividend / divisor   
           
            !amount of ice available to evaporate in time DT :
            evidt = devidt * dt                      !what we could get


            ! logic here is identical to rain. could get fancy and make subroutine
            ! but duplication of code is easier. God bless the screen editor.

            !case ice available to evaporate is lower than the ice
            if (evidt.le.qi (l) ) then                !enough ice to last DT                                    
               
               !case the ice available to evaporate is higher than the deficit, the air will saturate
               if (sd.le.evidt) then                  !enough time to saturate
                    qi (l) = qi (l) - sd                 !remove ice
                    qv (l) = qsat (l)                    !saturate
                    t  (l) = t (l) - sd * src            !cool the parcel
                    return                                  !done                                          
               else                                   !not enough time                                
                    sd = sd-evidt                        !use what there is
                    qv (l) = qv (l) + evidt              !add vapor
                    t  (l) =  t (l) - evidt * src        !lose heat
                    qi (l) = qi (l) - evidt              !lose ice                                    
               endif                                  !go on,unsatisfied                                    
           
           !case ice available to evaporate is higher than the ice, all ice will evaporate. !not enough ice to last DT  :
           else                                          

               !but there is enough to saturate, use it:
               if (sd.le.qi (l) ) then                    
                    qv (l) = qsat (l)                
                    qi (l) = qi   (l) - sd  
                    t  (l) =  t   (l) - sd * src  
                    return  
               else                                   !not enough to sat
                    sd = sd-qi (l)  
                    qv (l) = qv (l) + qi (l)  
                    t  (l) = t (l) - qi (l) * src             
                    qi (l) = 0.0                      !all gone
               endif                                  !on to better things
               
           endif    !finished with ice
         ENDIF !finished with the SD decision
          


return
end subroutine evaporate
!=================================================================================================================


!=================================================================================================================
subroutine sublimate

     !**********************************************************************************
     !vapor to ice: use equation OT22
     !**********************************************************************************


     !**********************************************************************************
     implicit none
     !**********************************************************************************

          real(kind=kind_phys), parameter :: &
               heatsubl = 2834., & ! from global parameters
               cp = 1.004          ! from global parameters
          real(kind=kind_phys), parameter :: &
               src = heatsubl / cpd !
          real(kind=kind_phys), parameter :: tfreeze = 269.3
          real(kind=kind_phys) :: &
               dtsubh,   & !
               dividend, & !
               divisor,  & !
               subl        !


          !call mpas_log_write('    --- enter subroutine sublimate:')

     !*****************************************************************
     ! calculations for sublimation rate
     !*****************************************************************

          dtsubh = 0.  
    
      !selection criteria for sublimation
          if (t (l)  .gt. tfreeze  ) return  !sublimation only occurs when  the temperature is lower then the freeze temperature and        
          if (qv (l) .le. qsat (l) ) return  !when there is excess of vapor
    
    
      !from (OT); correction factors for units applied
          dividend = ( (1.e6 / rho (l) ) **0.475) * (qv (l) / qsat (l) - 1.) * (qi (l) **0.525) * 1.13
          divisor = 7.e5 + 4.1e6 / (10. * est (l) )                                         
          dtsubh = ABS (dividend / divisor)   !sublimation rate
          subl = dtsubh * dt                  !and amount possible for a given dt


     !*****************************************************************
     !sublimation process
     !*****************************************************************

          !     again check the possibilities

          if (subl.lt.qv (l) ) then  
               qv (l) = qv (l) - subl             !lose vapor
               qi (l) = qi (l) + subl             !gain ice
               t  (l) = t  (l) + subl * src       !energy change, warms air
               
               return  
          else 
          !case when subl>= qv, all qv will sublimate                                     
               qi (l) = qv (l)                    !use what there is
               t  (l) = t (l) + qv (l) * src      !warm the air
               qv (l) = 0.0

          endif


 return
end subroutine sublimate
!=================================================================================================================


!=================================================================================================================
subroutine glaciate

     !***********************************************************************************
     !    Conversion of rain to ice
     !    uses equation OT 16, simplest. 
     !    correction from W not applied, but
     !    vapor pressure differences are supplied.
     !***********************************************************************************


     !**********************************************************************
     implicit none
     !**********************************************************************


          real(kind=kind_phys), parameter :: &
               heatfus = 334., & !
               cp = 1.004        ! from global parameters
               
          real(kind=kind_phys), parameter :: &
               frc = heatfus / cpd,  & ! 
               tfreeze =  269.3        !
          real(kind=kind_phys), parameter :: glconst = 0.025   !glaciation time constant, 1/sec
          real(kind=kind_phys) :: dfrzh


          !call mpas_log_write('    --- enter subroutine glaciate:')

     !*****************************************************************
     ! calculations for glaciation rate (rate of mass gain in ice)
     !*****************************************************************
      
      !rate of mass gain in ice
          dfrzh = 0.  

      !selection criteria for glaciation

          if (qh (l) .le. 0.    ) return     !case rain=0
          if (qv (l) .lt. qsat (l) ) return  !case                                       
          if (t  (l) .gt. tfreeze  ) return  
          !      NT=TMELT-T(L)          
          !      IF (NT.GT.50) NT=50 
              
      !from (OT 16): 
          dfrzh = dt * glconst * qh (l)    ! from OT(16)

     !*****************************************************************
     !glaciation process
     !*****************************************************************

          if (dfrzh.lt.qh (l) ) then 
               qi (l) = qi (l) + dfrzh  
               qh (l) = qh (l) - dfrzh  
               t  (l) = t  (l) + frc * dfrzh  !warms air 

               return  
          else
          !case dfrzh>= qh  
               qi (l) = qi (l) + qh (l)  
               t  (l) = t  (l) + frc * qh (l)  
               qh (l) = 0.0

          endif

 return
end subroutine glaciate
!=================================================================================================================

!=================================================================================================================
subroutine melt

     !makes water out of ice

     !**********************************************************************************
     implicit none
     !**********************************************************************************
                                           
          real(kind=kind_phys), parameter :: &
               frc = 332.27, & !
               tmelt = 273., & ! 
               f0 = 0.75       ! ice velocity factor
          real(kind=kind_phys) :: dtmelt     ! conversion,ice to rain

          !call mpas_log_write('    --- enter subroutine melt:')

     !*****************************************************************
     ! calculations for melting rate 
     !*****************************************************************
      
      !conversion,ice to rain
          dtmelt = 0.

      !selection criteria for melting
  
          if (qi (l) .le. 0.0  ) return  !case no ice
          if (t (l)  .lt. tmelt) return  !case temperature is not high enough
    
              
      !from (OT 23,24), after Mason,1956: 
         
          dtmelt = dt * (2.27 / rho (l) ) * cvi(l) * (t (l) - tmelt) * ( (rho(l) * qi (l) * 1.e-6) **0.525) * (f0** ( - 0.42) )
          


     !*****************************************************************
     !melting process
     !*****************************************************************

          if (dtmelt.lt.qi (l) ) then  
               qh (l) = qh (l) + dtmelt  
               qi (l) = qi (l) - dtmelt  
               t  (l) = t  (l) - frc * dtmelt     !cools air

               return  
          else  
               qh (l) = qh (l) + qi (l)   !get all there is to get
               t  (l) = t (l) - frc * qi (l)  
               qi (l) = 0.0  
               
          endif
   
 return
end subroutine melt
!=================================================================================================================

!=================================================================================================================
subroutine convert

     !*************************************************************************************************
     ! Accretion and Autoconversion
     ! Simulates the process of autoconversion (transformation of cloudwater to rainwater)
     ! and accretion (agregation of smaller cloud droplets to large ones by rainwater)
     !*************************************************************************************************

     !*************************************************************************************************
      implicit none
     !*************************************************************************************************

          real(kind=kind_phys),  parameter ::  ak1 = 0.001    !conversion rate constant
          real(kind=kind_phys),  parameter ::  ak2 = 0.0052   !collection (accretion) rate
          real(kind=kind_phys),  parameter ::  th  = 0.5      !Kessler threshold (minumum cloud water for autoconversion)
          integer,   parameter ::  iconv = 1      !(iconv =1 for Berry autoconversion, iconv =0 for Kessler autoconversion)                                     
    
          !real(kind=kind_phys), parameter :: ANBASE =  50.!*1.e+6 !Berry-number at cloud base #/m^3(maritime) !can I delete?
          real(kind=kind_phys),  parameter :: anbase =100000.!*1.e+6 !Berry-number at cloud base #/m^3(continental) (droplet concentration at cloud base)
    
          ! na formulacao abaixo use o valor em #/cm^3  
          !real(kind=kind_phys), parameter :: bdisp = 0.366       !Berry--size dispersion (maritime)
          real(kind=kind_phys),  parameter :: bdisp = 0.146       !Berry--size dispersion (continental). droplet size distribution 
          real(kind=kind_phys),  parameter :: tfreeze = 269.3  !ice formation temperature
          real(kind=kind_phys) :: &
               accrete, & !rate of accretion
               con,     & !rate of autoconversion
               q,       & !
               h,       & !
               total      ! total rate of accretion plus autoconversion


          !call mpas_log_write('    --- enter subroutine convert:')

     !*****************************************************************
     ! calculations for melting rate 
     !*****************************************************************
      
                
      !selection criteria for convert
          
          if (t (l)  .le. tfreeze) return              !process not allowed above ice
          if (qc (l) .eq. 0.    ) return               !process dont occurs if there's no cloud water  
          !srf IF (QC (L) .lt. 1.e-3     ) RETURN      

              
      !equations for accretion and autoconversion rates: 
         
          accrete = 0.
          con = 0. 
          q = rho (l) * qc (l)  !cloud water concentration
          h = rho (l) * qh (l)  !rain water concentration


          !-- accretion, Kessler (cloud droplets are collected by rain drops according to ak2)
          if (qh (l) .gt. 0.    ) accrete = ak2 * q * (h**.875)  


          !-- case Berry scheme for autoconversion
          if (iconv.ne.0) then  
               !old   bc1 = 120.  
               !old   bc2 = .0266 * anbase * 60.  
               !old   con = bdisp * q * q * q / (bc1 * q * bdisp + bc2)    
               con = q*q*q*bdisp/(60.*(5.*q*bdisp+0.0366*anbase))
          else 
          !-- case Kessler scheme for autoconversion     
               !   con = ak1 * (q - th)   !kessler autoconversion rate   
               !   IF (con.LT.0.0) con = 0.0   !havent reached threshold
               con = max(0.,ak1 * (q - th)) ! versao otimizada
          endif


          !--total rate
          total = (con + accrete) * dt / rho (l)  



     !*****************************************************************
     !conversion process
     !*****************************************************************

          if (total.lt.qc (l) ) then  
               qc (l) = qc (l) - total  
               qh (l) = qh (l) + total    !no phase change involved
               return 
          else       
               qh (l) = qh (l) + qc (l)    !uses all there is
               qc (l) = 0.0
          endif



 return
end subroutine convert
!=================================================================================================================


!=================================================================================================================
subroutine fallpart(m1)

     !************************************************************
     implicit none
     !************************************************************

          !     XNO=1.E7  [m**-4] median volume diameter raindrop,Kessler
          !     VC = 38.3/(XNO**.125), median volume fallspeed eqn., Kessler
          !     for ice, see (OT18), use F0=0.75 per argument there. rho*q
          !     values are in g/m**3, velocities in m/s

          integer, intent(in) :: m1 ! Index of upper boundary plume
          integer :: k ! loop counter model level
          real(kind=kind_phys) :: &
               vtc,  & ! NP: not used at all
               dfhz, & !
               dfiz, & !
               dz1     !
          
          real(kind=kind_phys), parameter :: &
               vconst = 5.107387, & !
               eps = 0.622,       & !
               f0 = 0.75            !  

          !call mpas_log_write('    --- enter subroutine fallpart:')

     !************************************************************
     !
     !************************************************************

          do k=2,m1-1
               !  vtc = vconst * rho (k) **.125   ! median volume fallspeed (KTable4) . this is not comment in the original. ?? why niko commented??                             
               !  hydrometeor assembly velocity calculations (k table4)
               !  vth(k)=-vtc*qh(k)**.125  !median volume fallspeed, water            
               vth (k) = - 4.      !small variation with qh
               vhrel = w (k) + vth (k)  !relative to surrounding cloud
               !  rain ventilation coefficient for evaporation
               cvh(k) = 1.6 + 0.57e-3 * (ABS (vhrel) ) **1.5
               !  vti(k)=-vtc*f0*qi(k)**.125    !median volume fallspeed,ice             
               vti (k) = - 3.                !small variation with qi
               virel = w (k) + vti (k)       !relative to surrounding cloud
               !  ice ventilation coefficient for sublimation
               cvi(k) = 1.6 + 0.57e-3 * (ABS (virel) ) **1.5 / f0  
          enddo ! 

          !NP: split loop for better vectorization. should i leave in this way? original its different


          do k=2,m1-1
               vhrel = w(k)  + vth (k)  !relative to surrounding cloud
               virel = w(k)  + vti(k)   !relative to surrounding cloud
               
               if (vhrel.ge.0.0) then  
                    dfhz=qh(k)*(rho(k  )*vth(k  )-rho(k-1)*vth(k-1))/rho(k-1)
                else  
                    dfhz=qh(k)*(rho(k+1)*vth(k+1)-rho(k  )*vth(k  ))/rho(k)
               endif

               if (virel.ge.0.0) then  
                    dfiz=qi(k)*(rho(k  )*vti(k  )-rho(k-1)*vti(k-1))/rho(k-1)
                else  
                    dfiz=qi(k)*(rho(k+1)*vti(k+1)-rho(k  )*vti(k  ))/rho(k)
               endif

               dz1=zm(k)-zm(k-1)

               qht(k) = qht(k) - dfhz / dz1 !hydrometeors don't

               qit(k) = qit(k) - dfiz / dz1  !nor does ice? hail, what about
          enddo

end subroutine fallpart
!=================================================================================================================


!=================================================================================================================
subroutine tend0_plumerise


     implicit none

          wt   (1:nm1)  = 0.
          tt   (1:nm1)  = 0.
          qvt  (1:nm1)  = 0.
          qct  (1:nm1)  = 0.
          qht  (1:nm1)  = 0.
          qit  (1:nm1)  = 0.
          vel_t(1:nm1)  = 0.
          rad_t(1:nm1)  = 0.
          sct  (1:nm1)  = 0.

     !call mpas_log_write('     --- enter subroutine tend0_plumerise:')

end subroutine tend0_plumerise
!=================================================================================================================


!=================================================================================================================
subroutine vel_advectc_plumerise(m1,wc,wt,rho,dzm)

     !**************************************************************
     implicit none
     !*************************************************************

          integer, intent(in) :: m1 ! model top index set in the loop of makeplume
          real(kind=kind_phys), dimension(m1), intent(in) :: &
               wc,  & ! set o zero in initial exept wc(1) = w(1)
               dzm, & ! inverse of full level vertical grid
               rho    ! air parcel dry air density g/m**3
          real(kind=kind_phys), dimension(m1), intent(inout) :: wt ! advection contribution - W tendency
          

          !-- local
          integer :: k ! loop counter
          real(kind=kind_phys), dimension(m1) :: flxw ! flux of w
          real(kind=kind_phys), dimension(m1) :: dn0 ! var local - converte de cgs para mks
          real(kind=kind_phys) :: c1z ! 

          !call mpas_log_write('    --- enter subroutine vel_advectc_plumerise') 

     !**************************************************************
     ! Flux of vertical velocity 
     !*************************************************************


      !-- converts units
          dn0(1:m1)=rho(1:m1)*1.e-3 !  converte de cgs(cm,g,sec) para mks (m,kg,sec). dn0=[kg air/m2.s]
                
      !-- computes the vertical velocity flux

          flxw(1) = wc(1) * dn0(1) 
      
          !For each point k>=2, it computes the flux as the product of wc(k) and the average density between levels k and k+1,
          !using a central difference method to account for density changes.
          do k = 2,m1-1
               flxw(k) = wc(k) * .5 * (dn0(k) + dn0(k+1)) !this gives me the flux in zmid or in in the interface??
          enddo


     !**************************************************************
     ! Computes advection contribution to W tendency 
     !*************************************************************

     
          c1z = .5

          do k = 2,m1-2
               wt(k) = wt(k)  &
                     + c1z * dzm(k) / (dn0(k) + dn0(k+1))             &
                     * ( (flxw(k)   + flxw(k-1))  * (wc(k) + wc(k-1)) &
                     -   (flxw(k)   + flxw(k+1))  * (wc(k) + wc(k+1)) &
                     +   (flxw(k+1) - flxw(k-1))  * 2. * wc(k)  )
          enddo

 return
end subroutine vel_advectc_plumerise
!=================================================================================================================



!=================================================================================================================
subroutine scl_advectc_plumerise (mzp)


     !*******************************************************************************
     implicit none
     !******************************************************************************

          integer, intent(in) :: mzp ! Index of upper plume boundary
          real(kind=kind_phys) :: dtlto2 ! 1/2 dt
          integer :: k ! loop counter model level

          !call mpas_log_write('    --- enter subroutine scl_advectc_plumerise')

     !*******************************************************************************
     ! Advect  scalars
     !******************************************************************************

          dtlto2   = .5 * dt
          vt3dc(1) = (w(1) + wc(1)) * dtlto2 * rho(1)*1.e-3 !converte de CGS p/ MKS
          vt3df(1) = .5 * (w(1) + wc(1)) * dtlto2 * dzm(1)

          do k = 2,mzp
               vt3dc(k) =  (w(k) + wc(k)) * dtlto2 *.5 * (rho(k) + rho(k+1))*1.e-3
               vt3df(k) =  (w(k) + wc(k)) * dtlto2 *.5 *  dzm(k)
          enddo

          do k = 1,mzp
               vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
               vctr2(k) = (zm(k)   - zt(k)) * dzm(k)
               vt3dk(k) = dzt(k) /(rho(k)*1.e-3)
          enddo

          !      scalarp => scalar_tab(n,ngrid)%var_p 
          !      scalart => scalar_tab(n,ngrid)%var_t 


     !**********************************************************************************
     ! Advection tendencies 
     !********************************************************************************** 
      
          !- Temperature advection tendency (tt)
               scr1=t
               call fa_zc_plumerise(mzp,t,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)
               call advtndc_plumerise(mzp,t,scr1(1),tt)

          !- water vapor advection tendency (qvt)
               scr1=qv
               call fa_zc_plumerise(mzp,qv,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk (1),vctr1,vctr2)
               call advtndc_plumerise(mzp,qv,scr1(1),qvt)



          !- liquid advection tendency (qct)
               scr1=qc
               call fa_zc_plumerise(mzp,qc,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)
               call advtndc_plumerise(mzp,qc,scr1(1),qct)


          !- ice advection tendency (qit)
               scr1=qi
               call fa_zc_plumerise(mzp,qi,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)
               call advtndc_plumerise(mzp,qi,scr1(1),qit)


          !- hail/rain advection tendency (qht)
               !   IF(ak1 > 0. .or. ak2 > 0.) THEN 
               scr1=qh
               call fa_zc_plumerise(mzp,qh,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)
               call advtndc_plumerise(mzp,qh,scr1(1),qht)
               !   ENDIF


          !- horizontal wind advection tendency (vel_t) !horizontal wind is considered scalar?? so o modulo
               scr1=vel_p
               call fa_zc_plumerise(mzp,vel_p,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)
               call advtndc_plumerise(mzp,vel_p,scr1(1),vel_t)


          !- vertical radius transport
               scr1=rad_p
               call fa_zc_plumerise(mzp,rad_p,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)
               call advtndc_plumerise(mzp,rad_p,scr1(1),rad_t)


          !- gas/particle advection tendency (SCT)  !I dont understand this part here??? comment
               scr1=SC
               call fa_zc_plumerise(mzp,sc,scr1(1),vt3dc(1),vt3df(1),vt3dg(1),vt3dk(1),vctr1,vctr2)
               call advtndc_plumerise(mzp,sc,scr1(1),sct)

 return
end subroutine scl_advectc_plumerise
!=================================================================================================================

!=================================================================================================================
subroutine fa_zc_plumerise(m1,scp,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2) !nao esta passando o argumento dt ????

     !****************************************************************************
     implicit none
     !****************************************************************************

          integer, intent(in) :: m1 ! Index upperboundary of plume
          
          real(kind=kind_phys), dimension(m1), intent(inout) :: &
               scr1, & !
               vt3dg   !
          real(kind=kind_phys), dimension(m1), intent(in) :: &
               vt3dk, & !
               vt3df, & !
               vt3dc, & !
               scp,   & !
               vctr1, & !
               vctr2    !

          !local
          integer :: k ! loop counter model level
          real(kind=kind_phys) :: dfact ! 1/2


          !call mpas_log_write('    --- enter subroutine fa_zc_plumerise')


     !*********************************************************************************
     ! Computes scalar flux 
     !*********************************************************************************

          dfact = .5

          ! Compute scalar flux VT3DG
          do k = 1,m1-1
               vt3dg(k) = vt3dc(k)                   &
                         * (vctr1(k) * scr1(k)        &
                         +  vctr2(k) * scr1(k+1)      &
                         +  vt3df(k) * (scr1(k) - scr1(k+1)))
          enddo

          ! Modify fluxes to retain positive-definiteness on scalar quantities.
          !    If a flux will remove 1/2 quantity during a timestep,
          !    reduce to first order flux. This will remain positive-definite
          !    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
          !    both fluxes are evacuating the box.

          do k = 1,m1-1
               if (vt3dc(k) .GT. 0.) THEN
                    if (vt3dg(k) * vt3dk(k) .GT. dfact * scr1(k)) THEN
                         vt3dg(k) = vt3dc(k) * scr1(k)
                    endif
               elseif (vt3dc(k) .LT. 0.) THEN
                    if (-vt3dg(k) * vt3dk(k+1) .GT. dfact * scr1(k+1)) THEN
                         vt3dg(k) = vt3dc(k) * scr1(k+1)
                    endif
               endif

          enddo

     !*********************************************************************************
     ! Computes flux divergence
     !*********************************************************************************

    
          do k = 2,m1-1
               scr1(k) = scr1(k)  &
                         + vt3dk(k) * ( vt3dg(k-1) - vt3dg(k) &
                         + scp  (k) * ( vt3dc(k)   - vt3dc(k-1)))
          enddo

 return
end subroutine fa_zc_plumerise
!=================================================================================================================

!=================================================================================================================
subroutine advtndc_plumerise (m1,scp,sca,sct) 
     
     !**********************************************************************************
     implicit none
     !********************************************************************************** 

          integer, intent(in) :: m1 ! Index upper boundary of plume
          real(kind=kind_phys), dimension(m1), intent(in) :: &
               scp, & !
               sca    !
          real(kind=kind_phys), dimension(m1), intent(inout) :: sct ! new tendency - NP: used for tendencies other then sct as well

          !local
          integer :: k ! loop counter model level
          real(kind=kind_phys) :: dtli ! 1/dt

          !call mpas_log_write('    --- enter subroutine advtndc_plumerise')


     !**********************************************************************************
     !
     !********************************************************************************** 
     
          dtli = 1. / dt
          
          do k = 2,m1-1
               sct(k) = sct(k) + (sca(k)-scp(k)) * dtli
          enddo

 return
end subroutine advtndc_plumerise
!=================================================================================================================


!=================================================================================================================
subroutine scl_advectc_plumerise2(mzp)

     !******************************************************************************
     implicit none
     !******************************************************************************

          integer,intent(in) :: mzp

          !- must be always real(8). 
          real(8), dimension(-5:nkp+5) :: u,v,p !jpereira20260430
          real(8), dimension(nkp) :: advtnd
          real(8) :: deltat,deltax
          integer i,k


          !====== for debug only
          ! w(:)=5.
          ! wc(:)=5.
          ! rho(:) =1.
          !rho_mks(:)=rho(:)!*0.001
          !v(:)=0.
          !====== 

          !call mpas_log_write('     -- enter subroutine scl_advectc_plumerise2')
  
     !**********************************************************************************
     ! Advection tendencies 
     !********************************************************************************** 

          do k=1,mzp
               v(k)=(w(k) + wc(k)) *.5
          enddo
  
          do k=-5,0
               v(k)=v(1)
          enddo
   
          deltax=dz  
          deltat=dt
  
          !-original definition of parameter 'p'
          !p=deltat/deltax
          p(1:nkp)=deltat*dzt(1:nkp) !  deltax=zm(K)-zm(k-1)
          p(-5:0) = deltat*dzt(1) !jpereira20260430

          
          !--- second term: 
               scr1=0.
               u(1:mzp)= rho(1:mzp)
               u(-5:0 )= rho(1)
               call bott(u,v,mzp,deltat,p,advtnd)
               scr1(1:mzp)=1.0*advtnd(1:mzp)

          
          !- temperature advection tendency (tt)
               u(1:mzp)= rho(1:mzp)*T(1:mzp)
               u(-5:0) =u(1)
               call bott(u,v,mzp,deltat,p,advtnd)
               tt(2:mzp) = tt(2:mzp)+advtnd(2:mzp)/rho(2:mzp) - &
                         (t (2:mzp)/rho(2:mzp))*scr1(2:mzp)


          !- water vapor advection tendency (qvt)
               u(1:mzp)= rho(1:mzp)*qv(1:mzp)
               u(-5:0) = u(1)
               call bott(u,v,mzp,deltat,p,advtnd)
               qvt(2:mzp) = qvt(2:mzp)+advtnd(2:mzp)/rho(2:mzp) - &
                              (qv (2:mzp)/rho(2:mzp))*scr1(2:mzp)


          !- horizontal wind advection tendency (vel_t)
               u(1:mzp)= rho(1:mzp)*vel_p(1:mzp)
               u(-5:0)=u(1)
               call bott(u,v,mzp,deltat,p,advtnd)
               vel_t(2:mzp) = vel_t(2:mzp)+advtnd(2:mzp)/rho(2:mzp) - &
                              ( vel_p(2:mzp)/rho(2:mzp))*scr1(2:mzp)


          !- vertical radius transport
               u(1:mzp)= rho(1:mzp)*rad_p(1:mzp)
               u(-5:0)=u(1)
               call bott(u,v,mzp,deltat,p,advtnd)
               rad_t(2:mzp) = rad_t(2:mzp)+advtnd(2:mzp)/rho(2:mzp) - &
                              ( rad_p(2:mzp)/rho(2:mzp))*scr1(2:mzp)



          if(do_microphysics) then  !should include this statement in the scl_advectc_plumerise??? and why the temperature is not include in the loop, since we have change in temperature from the microphysics too???

               !- liquid advection tendency (qct)
                    u(1:mzp)= rho(1:mzp)*QC(1:mzp)
                    u(-5:0)=u(1)
                    call bott(u,v,mzp,deltat,p,advtnd)
                    qct(2:mzp) = qct(2:mzp)+advtnd(2:mzp)/rho(2:mzp) - &
                                   (qc (2:mzp)/rho(2:mzp))*scr1(2:mzp)



               !- ice advection tendency (QIT)
                    u(1:mzp)= rho(1:mzp)*QI(1:mzp)
                    u(-5:0)=u(1)
                    call bott(u,v,mzp,deltat,p,advtnd)
                    qit(2:mzp) = qit(2:mzp)+advtnd(2:mzp)/rho(2:mzp) - &
                                   (qi (2:mzp)/rho(2:mzp))*scr1(2:mzp)


               
               !- hail advection tendency (QHT)
                    u(1:mzp) = rho(1:mzp)*QH(1:mzp)
                    u(-5:0 ) = u(1)
                    call bott(u,v,mzp,deltat,p,advtnd)
                    qht(2:mzp) = qht(2:mzp)+advtnd(2:mzp)/rho(2:mzp) - &
                                   (qh (2:mzp)/rho(2:mzp))*scr1(2:mzp)

          endif

     return 
          
          !- gas/particle advection tendency (SCT)
               u(1:mzp)= rho(1:mzp)*sc(1:mzp)
               u(-5:0)=u(1)
               call bott(u,v,mzp,deltat,p,advtnd)
               sct(2:mzp) = sct(2:mzp)+advtnd(2:mzp)/rho(2:mzp) - &
                              (sc (2:mzp)/rho(2:mzp))*scr1(2:mzp)

 return 
end subroutine scl_advectc_plumerise2

!===================================================================================================================


!===================================================================================================================
subroutine bott(u,v,k,deltat,p,advtnd)

     !****************************************************************
     !     Esta subrotina avalia a adveccao segundo o   esquema    
     !      de Bott(1989) - v. Chlond (1994) - "Locally Modified   
     !      Version of Bott's Advection Scheme" - MON. WEA. REV.   
     !      122,111-125                                            
     !
     !------ u: funcao a ser advectada ------------------------------
     !------ v: campo de velocidades, interpolado para os  pontos ---
     !------ intermediarios da grade - v(i) = 1/2 (vo(i)+vo(i+1)) ---
     !------ k: numero de passos da malha espacial ------------------
     !------ deltat: passo temporal ---------------------------------
     !------ p: deltat/passo espacial -------------------------------
     !------ ui, uf: valores de fronteira fixos ---------------------
     !------ lflag: 1 - algoritmo de 2a. ordem ----------------------
     !------        2 - algoritmo de 4a. ordem ----------------------
     !
     !****************************************************************


     !***************************************************************
     implicit none
     !***************************************************************
 
          ! real*8(a-h,o-z)
          integer, intent(in) :: k
          !real(kind=kind_phys),intent(in) :: deltat  !was defined as real(8) in the subroutine scl_advectc_plumerise2.
          real(8),intent(in) ::  deltat
          
          !-- must always be real*8  
          !real(8),intent(in) , dimension(k) ::  p !jpereira20260430
          real(8),intent(out), dimension(k) ::  advtnd
          real(8),dimension(-5:nkp+5)       :: u,v,cp,cm,aip,aim,beta,flux,p !jpereira20260430
          real(8),dimension(-5:nkp+5,0:8)   :: a
          real(8) :: eps, xzero,xum
          integer :: i,j,aj
 


          eps=1.e-20; xzero=0.; xum=1.
          !advtnd(:) = 0.
          !a(:,:)=0.
          !p(1:k)=deltat*dzt(1:k) ! dzt = 1/ (zm(k)-zm(k-1))
          !p(:)=deltat/deltax


          do i=0,k
               aip(i)=0.
               aim(i)=0.
          enddo

          !call mpas_log_write('    --- enter subroutine bott')




     !***************************************************************
     ! Calculo dos coeficientes
     !***************************************************************



          do i=0,k
               if (lflag.eq.0) then
                    a(i,0)=u(i)
               endif
               if (lflag.eq.1) then
                    a(i,0)=u(i)
                    a(i,1)=.5*(u(i+1)-u(i-1))
               endif
               if (lflag.eq.2) then
                    a(i,0)=-(u(i+1)-26.*u(i)+u(i-1))/24.
                    a(i,1)=.5*(u(i+1)-u(i-1))
                    a(i,2)=.5*(u(i+1)-2.*u(i)+u(i-1))
               endif
               if (lflag.eq.3) then
                    a(i,0)=-(u(i+1)-26.*u(i)+u(i-1))/24.
                    a(i,1)=(-5.*u(i+2)+34.*u(i+1)-34.*u(i-1)+5.*u(i-2))/48.
                    a(i,2)=.5*(u(i+1)-2.*u(i)+u(i-1))
                    a(i,3)=(u(i+2)-2.*u(i+1)+2.*u(i-1)-u(i-2))/12.
               endif
               if (lflag.eq.4) then
                    a(i,0)=(9.*u(i+2)-116.*u(i+1)+2134.*u(i)-116.*u(i-1)+9.*u(i-2))/1920.
                    a(i,1)=(-5.*u(i+2)+34.*u(i+1)-34.*u(i-1)+5.*u(i-2))/48.
                    a(i,2)=(-3.*u(i+2)+36.*u(i+1)-66.*u(i)+36.*u(i-1)-3.*u(i-2))/48.
                    a(i,3)=(u(i+2)-2.*u(i+1)+2.*u(i-1)-u(i-2))/12.
                    a(i,4)=(u(i+2)-4.*u(i+1)+6.*u(i)-4.*u(i-1)+u(i-2))/24.
               endif
               if (lflag.eq.5) then
                    a(i,0)=(9.*u(i+2)-116.*u(i+1)+2134.*u(i)-116.*u(i-1)+9.*u(i-2))/1920.
                    a(i,1)=(-259.*u(i-3)+2236.*u(i-2)-9455.*u(i-1)+9455.*&
                         u(i+1)-2236.*u(i+2)+259.*u(i+3))/11520.
                    a(i,2)=(-u(i+2)+12.*u(i+1)-22.*u(i)+12.*u(i-1)-u(i-2))/16.
                    a(i,3)=(7.*u(i-3)-52.*u(i-2)+83.*u(i-1)-83.*u(i+1)+52.&
                         *u(i+2)-7.*u(i+3))/288.
                    a(i,4)=(u(i+2)-4.*u(i+1)+6.*u(i)-4.*u(i-1)+u(i-2))/24.
                    a(i,5)=(-u(i-3)+4.*u(i-2)-5.*u(i-1)+5.*u(i+1)-4.*u(i+&
                         2)+u(i+3))/240.
               endif
               if (lflag.eq.6) then
                    a(i,0)=(-75.*u(i-3)+954.*u(i-2)-7621.*u(i-1)+121004.*&
                         u(i)-7621.*u(i+1)+954.*u(i+2)-75.*u(i+3))/107520.
                    a(i,1)=(-259.*u(i-3)+2236.*u(i-2)-9455.*u(i-1)+9455.*&
                         u(i+1)-2236.*u(i+2)+259.*u(i+3))/11520.
                    a(i,2)=(37.*u(i-3)-462.*u(i-2)+3435.*u(i-1)-6020.*u(i)&
                         +3435.*u(i+1)-462.*u(i+2)+37.*u(i+3))/3840.
                    a(i,3)=(7.*u(i-3)-52.*u(i-2)+83.*u(i-1)-83.*u(i+1)+52.&
                         *u(i+2)-7.*u(i+3))/288.
                    a(i,4)=(-5.*u(i-3)+54.*u(i-2)-171.*u(i-1)+244.*u(i)-&
                         171.*u(i+1)+54.*u(i+2)-5.*u(i+3))/576.
                    a(i,5)=(-u(i-3)+4.*u(i-2)-5.*u(i-1)+5.*u(i+1)-4.*u(i+&
                         2)+u(i+3))/240.
                    a(i,6)=(u(i-3)-6.*u(i-2)+15.*u(i-1)-20.*u(i)+15.*u(i+&
                         1)-6.*u(i+2)+u(i+3))/720.    
               endif
               if (lflag.eq.7) then
                    a(i,0)=(-75.*u(i-3)+954.*u(i-2)-7621.*u(i-1)+121004.*&
                              u(i)-7621.*u(i+1)+954.*u(i+2)-75.*u(i+3))/107520.
                    a(i,1)=(3229.*u(i-4)-33878.*u(i-3)+170422.*u(i-2)-&
                         574686.*u(i-1)+574686.*u(i+1)-170422.*u(i+2)+33878.*&
                              u(i+3)-3229.*u(i+4))/645120.
                    a(i,2)=(37.*u(i-3)-462.*u(i-2)+3435.*u(i-1)-6020.*u(i)&
                              +3435.*u(i+1)-462.*u(i+2)+37.*u(i+3))/3840.
                    a(i,3)=(-141.*u(i-4)+1406.*u(i-3)-6134.*u(i-2)+8614.*&
                              u(i-1)-8614.*u(i+1)+6134.*u(i+2)-1406.*u(i+3)+141.*&
                              u(i+4))/23040.
                    a(i,4)=(-5.*u(i-3)+54.*u(i-2)-171.*u(i-1)+244.*u(i)-&
                         171.*u(i+1)+54.*u(i+2)-5.*u(i+3))/576.
                    a(i,5)=(3.*u(i-4)-26.*u(i-3)+74.*u(i-2)-82.*u(i-1)+&
                         82.*u(i+1)-74.*u(i+2)+26.*u(i+3)-3.*u(i+4))/1920.
                    a(i,6)=(u(i-3)-6.*u(i-2)+15.*u(i-1)-20.*u(i)+15.*u(i+&
                         1)-6.*u(i+2)+u(i+3))/720.
                    a(i,7)=(-u(i-4)+6.*u(i-3)-14.*u(i-2)+14.*u(i-1)-14.*&
                         u(i+1)+14.*u(i+2)-6.*u(i+3)+u(i+4))/10080.
               endif
               if (lflag.eq.8) then
                    a(i,0)=(1225.*u(i-4)-17000.*u(i-3)+125884.*u(i-2)- &
                         800216.*u(i-1)+11702134.*u(i)-800216.*u(i+1)+125884.*&
                              u(i+2)-17000.*u(i+3)+1225.*u(i+4))/10321920.       
                    a(i,1)=(3229.*u(i-4)-33878.*u(i-3)+170422.*u(i-2)- &
                         574686.*u(i-1)+574686.*u(i+1)-170422.*u(i+2)+33878.*&
                              u(i+3)-3229.*u(i+4))/645120.                   
                    a(i,2)=(-3229.*u(i-4)+44480.*u(i-3)-323260.*u(i-2)+&
                         1912064.*u(i-1)-3260110.*u(i)+1912064.*u(i+1)-323260.*&
                              u(i+2)+44480.*u(i+3)-3229.*u(i+4))/1935360.  
                    a(i,3)=(-141.*u(i-4)+1406.*u(i-3)-6134.*u(i-2)+8614.*&
                              u(i-1)-8614.*u(i+1)+6134.*u(i+2)-1406.*u(i+3)+141.*&
                              u(i+4))/23040.
                    a(i,4)=(47.*u(i-4)-616.*u(i-3)+3908.*u(i-2)-10840.*&
                              u(i-1)+15002.*u(i)-10840.*u(i+1)+3908.*u(i+2)-616.*&
                              u(i+3)+47.*u(i+4))/27648.
                    a(i,5)=(3.*u(i-4)-26.*u(i-3)+74.*u(i-2)-82.*u(i-1)+&
                         82.*u(i+1)-74.*u(i+2)+26.*u(i+3)-3.*u(i+4))/1920.
                    a(i,6)=(-7.*u(i-4)+80.*u(i-3)-340.*u(i-2)+752.*u(i-1)-&
                              970.*u(i)+752.*u(i+1)-340.*u(i+2)+80.*u(i+3)-7.*&
                              u(i+4))/17280.
                    a(i,7)=(-u(i-4)+6.*u(i-3)-14.*u(i-2)+14.*u(i-1)-14.*&
                              u(i+1)+14.*u(i+2)-6.*u(i+3)+u(i+4))/10080.
                    a(i,8)=(u(i-4)-8.*u(i-3)+28.*u(i-2)-56.*u(i-1)+70.*&
                              u(i)-56.*u(i+1)+28.*u(i+2)-8.*u(i+3)+u(i+4))/40320.
               endif
          enddo


     !***************************************************************
     ! Calculo dos vetores (+) e (-) de numeros de Courant
     !***************************************************************

          do i=0,k
               cp(i)= dmax1(xzero,(p(i)*v(i)))
               cm(i)=-dmin1(xzero,(p(i)*v(i)))
          enddo


     !***************************************************************
     ! Calculo das integrais I(+,-)
     !***************************************************************


          do i=0,k
          do j=0,l
               aj=j
               aip(i)=aip(i)+a(i,j)*(1.-(1.-2.*cp(i))**(aj+1))/(aj+1.)/(2**(aj+1))
               aim(i)=aim(i)+a(i+1,j)*(-1)**j*(1.-(1.-2.*cm(i))**(aj+1))/(aj+1.)/(2**(aj+1))
          enddo
          enddo

          do i=0,k
               aip(i)=dmax1(aip(i),xzero)
               aim(i)=dmax1(aim(i),xzero)
          enddo

     !***************************************************************
     ! Calculo dos parametros beta
     !***************************************************************

          do i=0,k
               beta(i)=dmin1(xum,(u(i)/(dmax1((aip(i)+aim(i-1)),eps))))
          enddo

     !***************************************************************
     ! Calculo dos fluxos laterais
     !***************************************************************

          do i=0,k
               flux(i)=(beta(i)*aip(i)-beta(i+1)*aim(i))/p(i)
          enddo

     !***************************************************************
     ! Calculo da adveccao
     !***************************************************************


          do i=2,k
               !-- scalar updated
               !u(i)=u(i)-p*(flux(i)-flux(i-1))

               !srf- neste ponto tendencia de u =-p*(flux(i)-flux(i-1)) / deltat  
               !srf- isto e'  (du/dt)_adv=-p*(flux(i)-flux(i-1)) / deltat  
               advtnd(i) = -p(i)*(flux(i)-flux(i-1))/deltat  
          
          enddo

 return 
end  subroutine bott
!=================================================================================================================


!========================================================================================================
subroutine scl_misc(m1) 

     !*****************************************************************************
     implicit none
     !*****************************************************************************
          
          integer, intent(in) :: m1 ! index upper boundary of plume
          real(kind=kind_phys), parameter :: &
               g = 9.81, & ! from global parameters
               cp=1004.    ! from global parameters
          
          integer :: k ! loop counter model level
          real(kind=kind_phys) :: dmdtm ! (1/m)dm/dt

          !call mpas_log_write('    --- enter subroutine scl_misc')

     !*****************************************************************************
     ! Updates the tendencies due to entrainment ?
     !*****************************************************************************
      
          do k=2,m1-1
               
               wbar    = 0.5*(w(k)+w(k-1))  
                              
               !-- dry adiabat
               adiabat = - wbar * grav / cpd 

               !-- entrainment     
               dmdtm = 2. * alpha(k) * ABS (wbar) / radius (k)  != (1/m)dm/dt

               !-- tendency temperature = adv + adiab + entrainment
               tt(k) = tt(k) + adiabat - dmdtm * ( t  (k) -    te (k) ) 

               !-- tendency water vapor = adv  + entrainment
               qvt(k) = qvt(k)         - dmdtm * ( qv (k) - qvenv (k) )
               qct(k) = qct(k)         - dmdtm * ( qc (k)  )
               qht(k) = qht(k)         - dmdtm * ( qh (k)  )
               qit(k) = qit(k)         - dmdtm * ( qi (k)  )

               !-- tendency horizontal speed = adv  + entrainment
               vel_t(k) = vel_t(k)     - dmdtm * ( vel_p (k) - vel_e (k) )

               !-- tendency plume radius = adv  + entrainment
               rad_t(k) = rad_t(k)     + 0.5*dmdtm*(6./5.)*radius (k)

               !-- tendency gas/particle = adv  + entrainment
               sct(k) = sct(k)         - dmdtm * ( sc (k) -   sce (k) ) !this part is different from the original. why?

    enddo

end subroutine scl_misc
!======================================================================================================


!======================================================================================================
subroutine scl_dyn_entrain(m1)


     !**************************************************
     implicit none
     !**************************************************

          integer, intent(in) :: m1 ! index upper boundary plume = nm1
          integer :: k ! loop counter model level
          real(kind=kind_phys) :: dmdtm ! (1/m)dm/dt
          
          !call mpas_log_write('    --- enter subroutine scl_dyn_entrain')

     !*****************************************************************************
     ! Updates the tendencies due to dynamic entrainment ?
     !*****************************************************************************
     
          do k=2,m1-1
   
               !-- tendency horizontal radius from dyn entrainment
               rad_t(K) = rad_t(K)   + ABS((vel_e(k)-vel_p(k)))/pi

               !-- entrainment     
               dmdtm = (2./pi)  *  ABS(vel_e (k) - vel_p (k)) / radius (k)  
               
               !-- tendency horizontal speed  from dyn entrainment
               vel_t(k) = vel_t(k)     - dmdtm * ( vel_p (k) - vel_e (k) )

               !     if(vel_p (k) - vel_e (k) > 0.) cycle 

               !-- tendency temperature  from dyn entrainment
               tt(k) = tt(k)           - dmdtm * ( t (k) - te  (k) ) 

               !-- tendency water vapor  from dyn entrainment
               qvt(k) = qvt(k)         - dmdtm * ( qv (k) - qvenv (k) )
               qct(k) = qct(k)         - dmdtm * ( qc (k)  )
               qht(k) = qht(k)         - dmdtm * ( qh (k)  )
               qit(k) = qit(k)         - dmdtm * ( qi (k)  )

               !-- tendency gas/particle  from dyn entrainment
               sct(k) = sct(k)         - dmdtm * ( sc (k) - sce (k) ) !it is different from the original 

          enddo

end subroutine scl_dyn_entrain
!======================================================================================================

!======================================================================================================
subroutine damp_grav_wave(ifrom,nm1,deltak,zt,w,t,tt,qv,qh,qi,qc,te,pe,qvenv,vel_p,vel_t,vel_e) 


     !****************************************************************
     implicit none
     !****************************************************************

          integer, intent(in) :: &
               nm1,   & !
               ifrom, & !
               deltak   !

          real(kind=kind_phys), dimension(nm1), intent(inout):: &
               w,     & !
               t,     & !
               tt       !
          real(kind=kind_phys), dimension(nm1)::  dummy !
          real(kind=kind_phys), dimension(nm1), intent(in):: &
               qv,    & ! NP: not used here
               qh,    & ! NP: not used here
               qi,    & ! NP: not used here
               qc,    & ! NP: not used here
               te,    & !
               pe,    & ! NP: not used here
               qvenv, & ! NP: not used here
               zt,    & !
               vel_p, & ! NP: not used here
               vel_t, & ! NP: not used here
               vel_e    ! NP: not used here

          !call mpas_log_write('    --- enter subroutine damp_grav_wave')

     !************************************************************
     !
     !************************************************************

          if(ifrom==1) then
               call friction(ifrom,nm1,deltak,zt,t,tt,te) !again theres no dt or zm.?
               ! call friction(ifrom,nm1,deltak,dt,zt,zm,vel_p,vel_t,vel_e)
               ! call friction(ifrom,nm1,dt,zt,zm,qv,qvt,qvenv)
               return
          endif

          dummy(:) = 0.
          if (ifrom==2) call friction(ifrom,nm1,deltak,zt,w,dummy ,dummy) !again theres no dt or zm.?
          !call friction(ifrom,nm1,dt,zt,zm,qi,qit ,dummy)
          !call friction(ifrom,nm1,dt,zt,zm,qh,qht ,dummy)
          !call friction(ifrom,nm1,dt,zt,zm,qc,qct ,dummy)
 return
end subroutine damp_grav_wave
!======================================================================================================

!======================================================================================================
subroutine friction(ifrom,nm1,deltak,zt,var1,vart,var2)
    
     !***********************************************
     implicit none
     !************************************************

    integer, intent(in) :: &
         nm1,   & ! index upper boundary model level
         ifrom, & !
         deltak   !
    integer :: &
         k,   & ! loop counter model level
         nfpt, & ! NP: not used at all
         kf     ! Index
    real(kind=kind_phys), dimension(nm1), intent(inout) :: & 
         var1, & !
         vart    !
    real(kind=kind_phys), dimension(nm1), intent(in) :: & 
         var2, & !
         zt      !
    real(kind=kind_phys) :: &
         zmkf,   & !
         ztop,   & !
         distim, & !
         c1,     & !
         c2        !

     !call mpas_log_write('     --- enter subroutine friction')

     !*************************************************************
     !
     !************************************************************

          kf = nm1 - int(deltak)
          zmkf = zm(kf) 
          ztop = zm(nm1)

          distim = 60. 
          ! distim = MIN(3.*dt,60.) 


          c1 = 1. / (distim * (ztop - zmkf))
          c2 = dt * c1

          if(ifrom == 1) then  
               do k = nm1,2,-1
                    if (zt(k) .le. zmkf) cycle  !
                    vart(k) = vart(k)   + c1 * (zt(k) - zmkf)*(var2(k) - var1(k))
               enddo
          elseif(ifrom == 2) then
               do k = nm1,2,-1
                    if (zt(k) .le. zmkf) cycle  ! 
                    var1(k) =  var1(k) + c2 * (zt(k) - zmkf)*(var2(k) - var1(k))
               enddo
          endif

 return
end subroutine friction
!======================================================================================================

!======================================================================================================
subroutine visc_W(m1,deltak,kmt)


     !************************************************************
     implicit none
     !************************************************************

          integer, intent(in) :: &
               m1,     & !
               deltak, & ! NP: not used here
               kmt       !
          integer :: &
               k,  & !
               m2    ! NP: not used at all 
          real(kind=kind_phys) :: &
               dz1t,      & !
               dz1m,      & !
               dz2t,      & !
               dz2m,      & !
               d2wdz,     & !
               d2tdz,     & !
               d2qvdz,    & !
               d2qhdz,    & !
               d2qcdz,    & !
               d2qidz,    & !
               d2scdz,    & !
               d2vel_pdz, & !
               d2rad_dz     !
          integer :: &
               mI,   & !
               mF,   & !
               deltaM  !

          !call mpas_log_write('    --- enter subroutine visc_W')

     !********************************************************
     !
     !********************************************************

          mI = 2
          mF = min(m1,kmt-1)
          deltaM = 1

          do k=mI,mF,deltaM !v2
               !do k=2,m2-1 !orig
               dz1t   = 0.5*(zt(k+1)-zt(k-1))
               dz2t   = visc (k) / (dz1t * dz1t)  
               dz1m   = 0.5*(zm(k+1)-zm(k-1))
               dz2m   = VISC (k) / (DZ1M * dz1m)  
               d2wdz  = (w  (k + 1) - 2 * w  (k) + w  (k - 1) ) * dz2m  
               d2tdz  = (t  (k + 1) - 2 * t  (k) + t  (k - 1) ) * dz2t  
               d2qvdz = (qv (k + 1) - 2 * qv (k) + qv (k - 1) ) * dz2t  
               d2qhdz = (qh (k + 1) - 2 * qh (k) + qh (k - 1) ) * dz2t 
               d2qcdz = (qc (k + 1) - 2 * qc (k) + qc (k - 1) ) * dz2t  
               d2qidz = (qi (k + 1) - 2 * qi (k) + qi (k - 1) ) * dz2t  
               d2scdz = (sc (k + 1) - 2 * sc (k) + sc (k - 1) ) * dz2t 
               d2vel_pdz=(vel_p  (k + 1) - 2 * vel_p  (k) + vel_p  (k - 1) ) * dz2t
               d2rad_dz =(rad_p  (k + 1) - 2 * rad_p  (k) + rad_p  (k - 1) ) * dz2t

               wt(k) =   wt(k) + d2wdz
               tt(k) =   tt(k) + d2tdz       

               qvt(k) =  qvt(k) + d2qvdz 
               qct(k) =  qct(k) + d2qcdz
               qht(k) =  qht(k) + d2qhdz 
               qit(k) =  qit(k) + d2qidz     

               vel_t(k) =   vel_t(k) + d2vel_pdz
               rad_t(k) =   rad_t(k) + d2rad_dz

               sct(k) =  sct(k) + d2scdz !the original not calculates that. why??
              
          enddo

end subroutine visc_W
!======================================================================================================


!======================================================================================================
subroutine update_plumerise (m1,varn)

     !************************************************************
     implicit none
     !************************************************************

          integer, intent(in) :: m1 ! index upper boundary plume
          integer :: k ! loop counter model level
          character(len=*), intent(in) :: varn

          !call mpas_log_write('     --- enter subroutine update_plumerise')


     !************************************************************
     !
     !************************************************************

          if(varn == 'W') then
               do k=2,m1-1
                    w(k) =  w(k) +  wt(k) * dt  
               enddo
               return
          else 
               do k=2,m1-1
                    t(k)  =  t(k) +  tt(k) * dt  
                    qv(k) = qv(k) + qvt(k) * dt  
                    qc(k) = qc(k) + qct(k) * dt 
                    qh(k) = qh(k) + qht(k) * dt  
                    qi(k) = qi(k) + qit(k) * dt 
                    qv(k) = MAX(0., qv(k))   
                    qc(k) = MAX(0., qc(k))
                    qh(k) = MAX(0., qh(k))
                    qi(k) = MAX(0., qi(k))
                    vel_p(k) =  vel_p(k) + vel_t(k) * dt  
                    rad_p(k) =  rad_p(k) + rad_t(k) * dt  
                    sc(k)    =  sc(k)    + sct(k)  * dt !the original does not calculates this. ??
               enddo
          endif

end subroutine update_plumerise
!======================================================================================================


!======================================================================================================
subroutine hadvance_plumerise (iac,m1,wc,wt,w,mintime)

     !*********************************************
     implicit none
     !*********************************************

          integer :: k ! NP: not used here 
          integer, intent(in) :: &
               iac,  & ! case choice
               m1,   & ! index upper boudary plume
               mintime !
          real(kind=kind_phys), dimension(m1) :: dummy
          real(kind=kind_phys), dimension(m1), intent(inout) ::  &
               wc, & !
               wt, & !
               w     !
          real(kind=kind_phys) :: eps

          !call mpas_log_write('     --- enter subroutine hadvance_plumerise')

     !********************************************************************************
     !
     !********************************************************************************


          !     It is here that the Asselin filter is applied.  For the velocities
          !     and pressure, this must be done in two stages, the first when
          !     IAC=1 and the second when IAC=2.

          !Asselin??

          eps = .2
          if(mintime == 1) eps=0.5

          
          !For both IAC=1 and IAC=2, call PREDICT for U, V, W, and P.
          call predict_plumerise(m1,wc,w,wt,dummy,iac,2.*dt,eps)

 return
end subroutine hadvance_plumerise
!======================================================================================================


!======================================================================================================
subroutine predict_plumerise(npts,ac,ap,fa,af,iac,dtlp,epsu)

     !*********************************************************
     implicit none
     !*********************************************************

          integer, intent(in) :: &
               npts, & ! index upper boundary plume
               iac     !
          integer :: m ! loop counter model levels
          real(kind=kind_phys), intent(in) :: epsu
          real(kind=kind_phys), intent(in) :: dtlp ! not used at all
          real(kind=kind_phys), dimension(*), intent(inout) :: &
               ac, & !
               ap, & !
               fa, & ! not used at all
               af    !

          !call mpas_log_write('    --- enter subroutine predict_plumerise')


     !*************************************************************************
     !     For IAC=3, this routine moves the arrays AC and AP forward by
     !     1 time level by adding in the prescribed tendency. It also
     !     applies the Asselin filter given by:

     !              {AC} = AC + EPS * (AP - 2 * AC + AF)

     !     where AP,AC,AF are the past, current and future time levels of A.
     !     All IAC=1 does is to perform the {AC} calculation without the AF
     !     term present.  IAC=2 completes the calculation of {AC} by adding
     !     the AF term only, and advances AC by filling it with input AP
     !     values which were already updated in ACOUSTC.
     !************************************************************************

     !***********************************************************************
     !
     !***********************************************************************

          if (iac .eq. 1) then
               do m = 1,npts
                    ac(m) = ac(m) + epsu * (ap(m) - 2. * ac(m))
               enddo
               return
          elseif (iac .eq. 2) then
               do m = 1,npts
                    af(m) = ap(m)
                    ap(m) = ac(m) + epsu * af(m)
               enddo
               !elseif (iac .eq. 3) then
               !   do m = 1,npts
               !      af(m) = ap(m) + dtlp * fa(m)
               !   enddo
               !   if (ngrid .eq. 1 .and. ipara .eq. 0) call cyclic(nzp,nxp,nyp,af,'T')
               !   do m = 1,npts
               !      ap(m) = ac(m) + epsu * (ap(m) - 2. * ac(m) + af(m))
               !   enddo
          endif

          do m = 1,npts
               ac(m) = af(m)
          enddo

 return
end subroutine predict_plumerise
!=================================================================================================================


!=================================================================================================================
subroutine buoyancy_plumerise (m1, t, te, qv, qvenv, qh, qi, qc, wt, scr1)


     !*******************************************************************
     implicit none
     !*******************************************************************

          integer :: k ! loop counter
          integer, intent(IN) :: m1 ! index upper boundary plume
          real(kind=kind_phys), parameter :: &
               g = 9.8,     & ! from global parameters
               eps = 0.622, & !
               gama = 0.5     ! mass virtual coeff.
          real(kind=kind_phys), dimension(m1), intent(in) :: &
               t,     & !
               te,    & !
               qv,    & !
               qvenv, & !
               qh,    & !
               qi,    & !
               qc       !
          real(kind=kind_phys), dimension(m1), intent(inout) :: &
               wt,    & !
               scr1     !
          real(kind=kind_phys) :: &
               tv,     & !
               tve,    & !
               qwtotl, & !
               umgamai   !
          real(kind=kind_phys), parameter :: mu = 0.15 


          !call mpas_log_write('   --- enter subroutine buoyancy_plumerise')


     !***********************************************************************************
     !
     !***********************************************************************************

          !- orig
          umgamai = 1./(1.+gama) ! compensa a falta do termo de aceleracao associado `as
                                 ! das pertubacoes nao-hidrostaticas no campo de pressao

          !- new                 ! Siesbema et al, 2004
          !umgamai = 1./(1.-2.*mu)

          do k = 2,m1-1
               tv =   t(k) * (1. + (qv(k)   /eps))/(1. + qv(k)   )  !blob virtual temp.                                               
               tve = te(k) * (1. + (qvenv(k)/eps))/(1. + qvenv(k))  !and environment
               qwtotl = qh(k) + qi(k) + qc(k)                       ! qwtotl*g is drag
               scr1(k)= grav*  umgamai*( (tv - tve) / tve   - qwtotl) 
          enddo

          do k = 2,m1-2
               !srf- just for output
               rbuoy(k)=0.5*(scr1(k)+scr1(k+1)) !if not used, delete.

               wt(k) = wt(k)+0.5*(scr1(k)+scr1(k+1))
          enddo
          
          !srf- just for output
          rbuoy(1)=rbuoy(2) !if not used, delete.


end subroutine buoyancy_plumerise
!=================================================================================================================


!=================================================================================================================
subroutine entrainment(m1,w,wt,radius,alpha)


     !*******************************************************************
     implicit none
     !*******************************************************************

          integer :: k ! loop counter model level
          integer, intent(in) :: m1 ! index upper boundary plume
          real(kind=kind_phys), dimension(m1), intent(in) :: &
               w,      & !
               radius, & !
               alpha       
          real(kind=kind_phys), dimension(m1), intent(inout) :: wt ! tendency of w
          real(kind=kind_phys) :: &
               dmdtm,      & !
               wbar,       & !
               radius_bar, & !
               umgamai,    & !
               dyn_entr      ! diynamic entrainment
          real(kind=kind_phys), parameter :: &
               mu = 0.15, & !
               gama = 0.5   ! mass virtual coeff.


          !call mpas_log_write('    --- enter subroutine entrainment')

     !*******************************************************************
     !
     !*******************************************************************

          !- new - Siesbema et al, 2004
          !umgamai = 1./(1.-2.*mu)

          !- orig
          umgamai = 1./(1.+gama) ! compensa a falta do termo de aceleracao associado `as
                                       !  pertubacoes nao-hidrostaticas no campo de pressao

          !-- ALPHA/RADIUS(L) = (1/M)DM/DZ  (W 14a)
          do k=2,m1-1
               !-- for W: WBAR is only W(k)       
               wbar=w(k)          
               radius_bar = 0.5*(radius(k) + radius(k-1))
               dmdtm = umgamai * 2. * alpha(k) * ABS (wbar) / radius_bar  != (1/m)dm/dt

               !--  dmdtm*w(l) entrainment,
               wt(k) = wt(k)  - dmdtm*ABS (wbar) !- dmdtm*abs (wbar)*1.875*0.5 !different from original that dont use ABS here. but brams use.  ??

               !if(vel_p (k) - vel_e (k) > 0.) cycle

               !-   dynamic entrainment
               dyn_entr =  (2./pi)*0.5*ABS (vel_p(k)-vel_e(k)+vel_p(k-1)-vel_e(k-1)) /radius_bar
               wt(k) = wt(k)  - dyn_entr*ABS (wbar) !different from original that dont use ABS here. but brams use. ??

               !- entrainment acceleration for output only
                dwdt_entr(k) =  - dmdtm*ABS (wbar)- dyn_entr*ABS (wbar) !delete if not used. 
          ENDDO
    
          !- entraiment acceleration for output only
          dwdt_entr(1) = dwdt_entr(2)

end subroutine entrainment
!=================================================================================================================

!=================================================================================================================
subroutine set_flam_vert(zzcon,ztopmax,k1,k2,zbot_plume,ztop_plume) !ask saulo to see if this subroutine is correct. ??
  
     !**************************************************************
     implicit none
     !**************************************************************
          real(kind=kind_phys), intent(in)  :: zzcon(:)
          real(kind=kind_phys),    dimension(2), intent(in)  :: ztopmax
          integer, intent(out) ::&
               k1,  &  
               k2
          
          real(kind=kind_phys),  intent(out) :: zbot_plume, ztop_plume
          
          !local:
          integer :: k
          integer, dimension(2)  :: k_lim    ! store the last valid index for each imm

          !call mpas_log_write('    --- enter subroutine set_flam_vert')

     !**************************************************************
     ! Define MPAS dynamical level index for the plume
     !**************************************************************

          !This loop checks each element of zzcon at level k against the ztopmax(imm), which means
          !that checks the height of the MPAS dynamical levels against the plume injection 
          !height calculated at PRM dynamical levels. When the MPAS height is higher than the plume height,
          !the loop stop and save the value of the vertical index of MPAS dynamical level in k_lim.

          
          do imm=1,2
                              
               do k=1,nkp-1                    
                    if(zzcon(k) > ztopmax(imm) ) exit
               enddo
               
               k_lim(imm) = k
          enddo  


          k1=max(3,k_lim(1))
          k2=max(3,k_lim(2))
          
          if(k2 < k1) then
               !print*,'1: ztopmax k=',ztopmax(1), k1
               !print*,'2: ztopmax k=',ztopmax(2), k2
               k2 = k1
          endif

          !if need the plume height in m MSL need to some with the surface height
          zbot_plume = zzcon(k1)
          ztop_plume = zzcon (k2)

    
end subroutine set_flam_vert
!=================================================================================================================



!=================================================================================================================
real function esat(t)

     !****************
      implicit none
     !****************
      
      real(kind=kind_phys) :: t
      real(kind=kind_phys), parameter :: abz=273.15
      real(kind=kind_phys) :: tc

     !****************
     ! esat function
     !****************
      
      !esat(millibars),t(kelvin)
      tc=t-abz
      esat=6.1078*exp((17.2693882*tc)/(tc+237.3))

 return
end function  esat
!=================================================================================================================



!==================================================================================================================
 end module mpas_chemistry_plumerisemodel_run
!==================================================================================================================
