! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
! Update the aerosol concentrations on the MPAS grid  due to emissions distributed at PRM injH
! Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-05-20.
!
!================================================================================================================== 
 module mpas_chemistry_gocart2G_plumerisemodel
 
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types 
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_array, mpas_pool_get_dimension, mpas_pool_get_config

 use mpas_chemistry_gocart2G_vars,only: mpas_gocart2G, mpas_emis_gocart2G
 use mpas_chemistry_plumerisemodel_vars,only: mpas_prm

 use CA2G_bc_instance, only: bc_fhydrophobic => hydrophobic_fraction !0.8 
 use CA2G_br_instance, only: br_fhydrophobic => hydrophobic_fraction !0.5
 use CA2G_oc_instance, only: oc_fhydrophobic => hydrophobic_fraction !0.5
 

 
 implicit none
    
    private
    public:: gocart2G_plumerisemodel

 contains

!==================================================================================================================
 subroutine gocart2G_plumerisemodel(configs, state, time_lev)
!==================================================================================================================

   !*******************************************************************************
   !implicit none
   !*******************************************************************************
    
      !input arguments:
         type(mpas_pool_type),intent(in):: configs
         type(mpas_pool_type),intent(in):: state
         integer,intent(in):: time_lev
      
      
      !local variables and pointers:

         integer:: i,its,ite,j,jts,jte,k,kts,kte
         real(kind=RKIND),pointer:: dt
         

         !local arrays for mixing ratio due to emissions
         real(kind=RKIND),dimension(:,:,:),allocatable:: f_smold          ! Factor applied to smoldering phase    
         real(kind=RKIND),dimension(:,:,:),allocatable:: f_flame          ! Factor applied to flaming phase phase
         real(kind=RKIND),dimension(:,:,:),allocatable:: qbcphobic_srcprm ! BC hydrophobic mixing ratio due to BB emission source (kg.kg-1)
         real(kind=RKIND),dimension(:,:,:),allocatable:: qbcphilic_srcprm ! BC hydrophilic mixing ratio due to BB emission source (kg.kg-1)
         real(kind=RKIND),dimension(:,:,:),allocatable:: qbrphobic_srcprm ! BR hydrophobic mixing ratio due to BB emission source (kg.kg-1)
         real(kind=RKIND),dimension(:,:,:),allocatable:: qbrphilic_srcprm ! BR hydrophilic mixing ratio due to BB emission source (kg.kg-1)
         real(kind=RKIND),dimension(:,:,:),allocatable:: qocphobic_srcprm ! OC hydrophobic mixing ratio due to BB emission source (kg.kg-1)
         real(kind=RKIND),dimension(:,:,:),allocatable:: qocphilic_srcprm ! OC hydrophilic mixing ratio due to BB emission source (kg.kg-1)
         real(kind=RKIND),dimension(:,:,:),allocatable:: qso2_srcprm      ! SO2 mixing ratio due to BB emission source (kg.kg-1)
         real(kind=RKIND),dimension(:,:,:),allocatable:: qnh3_srcprm      ! NH3 mixing ratio due to BB emission source (kg.kg-1)
         real(kind=RKIND),dimension(:,:,:),allocatable:: qsoapbb_srcprm   ! SOA precursor mixing ratio due to BB emission source (kg.kg-1)


         !local arrays to received MPAS-GOCART2 species mixing ratio
         real(kind=RKIND),dimension(:,:,:),pointer:: scalars
         integer,pointer:: index_qbcphobic   !BCPHOBIC - black carbon hydrophobic mass mixing ratio (kg kg^{-1})
         integer,pointer:: index_qbcphilic   !BCPHILIC - black carbon hydrophilic mass mixing ratio (kg kg^{-1})
         integer,pointer:: index_qbrphobic   !BRPHOBIC - brown carbon hydrophobic mass mixing ratio (kg kg^{-1})
         integer,pointer:: index_qbrphilic   !BRPHILIC - brown carbon hydrophilic mass mixing ratio (kg kg^{-1})
         integer,pointer:: index_qocphobic   !OCPHOBIC - organic carbon hydrophobic mass mixing ratio (kg kg^{-1})
         integer,pointer:: index_qocphilic   !OCPHILIC - organic carbon hydrophilic mass mixing ratio (kg kg^{-1})
         integer,pointer:: index_qso2        !SO2 - sulfur dioxide mass mixing ratio (kg kg^{-1})
         integer,pointer:: index_qnh3        !NH3 - ammonia mass mixing ratio (kg kg^{-1})
         integer,pointer:: index_qsoapbb     !SOAPbb - lumped SOA precursor mass mixing ratio from biomass burning emissions (kg kg^{-1})

      
      !variables accessed directly trough types:
         !mpas_prm%flame_fraction         ! Fraction emitted during the flaming phase according to the vegetation category
         !mpas_prm%rho_dry                ! Dry air density (kg m^{-3})
         !mpas_prm%dz                     ! Layer thickness  (m)
         !mpas_emis_gocart2G%bc_biomass   ! Black carbon biomass burning emission flux (kg/m^2/s)
         !mpas_emis_gocart2G%br_biomass   ! Brown carbon biomass burning emission flux (kg/m^2/s)
         !mpas_emis_gocart2G%oc_biomass   ! Organic carbon biomass burning emission flux (kg/m^2/s)
         !mpas_emis_gocart2G%su_biomass   ! Sulfur Dioxide (SO2) biomass burning emission flux (kg/m^2/s)
         !mpas_emis_gocart2G%nh3_bb       ! Amonia (NH3) biomass burning emission flux (kg/m^2/s)        
         !mpas_emis_gocart2G%soap_biomass ! SOA Precursor biomass burning emission flux (kg/m^2/s) 
         !bc_fhydrophobic                 ! Initially BC hydrophobic portion = 0.8
         !br_fhydrophobic                 ! Initially BR hydrophobic portion = 0.5
         !oc_fhydrophobic                 ! Initially OC hydrophobic portion = 0.5
         !mpas_prm%k1_zbot                ! Index level on MPAS grid of lower plume injection height 
         !mpas_prm%k2_ztop                ! Index level on MPAS grid of upper plume injection height


   call mpas_log_write('--- enter subroutine gocart2G_plumerisemodel:')
   call mpas_log_write('    Update the bb aerosol mix ratio due to emissions at PRM injH')


   !***********************************************************************************
   ! Calculates Aerosol Mixing Ratio due to emissions distributed on PRM injH
   !***********************************************************************************
   
      !-------------------------------------------------------------------------------!
      !                                                                               !
      !                                             (emis flux)                       !
      !aerosol mix ratio = (flame_fraction) *  _______________________  * time step   !
      !                                           dry air density * dz                !
      !-------------------------------------------------------------------------------!
      

      !--------------------------------------------------------------------------------
      !Get dimensions from MPAS needed to run GOCART2G
      !--------------------------------------------------------------------------------
         
         its    = mpas_gocart2G%its
         ite    = mpas_gocart2G%ite
         jts    = mpas_gocart2G%jts
         jte    = mpas_gocart2G%jte
         kts    = mpas_gocart2G%kts
         kte    = mpas_gocart2G%kte
            
         call mpas_log_write('ITS = $i   ITE = $i',intArgs=(/its,ite/))
         call mpas_log_write('JTS = $i   JTE = $i',intArgs=(/jts,jte/))
         call mpas_log_write('KTS = $i   KTE = $i',intArgs=(/kts,kte/))
   
      !-------------------------------------------------------------------------------------------------
      !Allocates local arrays for receiving Mixing Ratio Tendency due to emissions
      !-------------------------------------------------------------------------------------------------
         
         if(.not.allocated(f_smold)  )    allocate(f_smold(its:ite,kts:kte,jts:jte) )  
         if(.not.allocated(f_flame)  )    allocate(f_flame(its:ite,kts:kte,jts:jte) )

         if(.not.allocated(qbcphobic_srcprm)  ) allocate(qbcphobic_srcprm(its:ite,kts:kte,jts:jte) )
         if(.not.allocated(qbcphilic_srcprm)  ) allocate(qbcphilic_srcprm(its:ite,kts:kte,jts:jte) )
         if(.not.allocated(qbrphobic_srcprm)  ) allocate(qbrphobic_srcprm(its:ite,kts:kte,jts:jte) )
         if(.not.allocated(qbrphilic_srcprm)  ) allocate(qbrphilic_srcprm(its:ite,kts:kte,jts:jte) )
         if(.not.allocated(qocphobic_srcprm)  ) allocate(qocphobic_srcprm(its:ite,kts:kte,jts:jte) )
         if(.not.allocated(qocphilic_srcprm)  ) allocate(qocphilic_srcprm(its:ite,kts:kte,jts:jte) )
         if(.not.allocated(qso2_srcprm)  )      allocate(qso2_srcprm     (its:ite,kts:kte,jts:jte) )
         if(.not.allocated(qnh3_srcprm)  )      allocate(qnh3_srcprm     (its:ite,kts:kte,jts:jte) )
         if(.not.allocated(qsoapbb_srcprm)  )   allocate(qsoapbb_srcprm  (its:ite,kts:kte,jts:jte) )

      
      !----------------------------------------------------------------------------------
      !Initialize Arrays for Mixing Ratio due to emissions
      !----------------------------------------------------------------------------------
      
         do j = jts, jte
            do i = its, ite
               do k = kts, kte
                  
                  f_smold(i,k,j)      = 0._RKIND
                  f_flame(i,k,j)      = 0._RKIND

                  qbcphobic_srcprm(i,k,j)  = 0._RKIND
                  qbcphilic_srcprm(i,k,j)  = 0._RKIND
                  qbrphobic_srcprm(i,k,j)  = 0._RKIND
                  qbrphilic_srcprm(i,k,j)  = 0._RKIND  
                  qocphobic_srcprm(i,k,j)  = 0._RKIND
                  qocphilic_srcprm(i,k,j)  = 0._RKIND
                  qso2_srcprm(i,k,j)       = 0._RKIND
                  qnh3_srcprm(i,k,j)       = 0._RKIND
                  qsoapbb_srcprm(i,k,j)    = 0._RKIND

               enddo
            enddo
         enddo      
      
      
      !----------------------------------------------------------------------------------
      ! Calculates Mixing Ratio due to Emission Source
      ! Reference: CAEmission routine in chemistry/gocart2G/GOCART2G_Process.F90
      !----------------------------------------------------------------------------------

         !--- Get Time Step
         call mpas_pool_get_config(configs,'config_dt',dt)
            
         !--- Calculates Mixing Ratio due to emission source   
         do j = jts, jte
            do i = its, ite
                        
               !Mixing ratio at surface (smoldering phase):  

               f_smold(i,kts,j) = (1 - mpas_prm%flame_fraction(i,j)) / (mpas_prm%rho_dry(i,kts,j) * mpas_prm%dz(i,kts,j))
               
               qbcphobic_srcprm(i,kts,j) = f_smold(i,kts,j) *  mpas_emis_gocart2G%bc_biomass(i,j)  *  bc_fhydrophobic      * dt 
               qbcphilic_srcprm(i,kts,j) = f_smold(i,kts,j) *  mpas_emis_gocart2G%bc_biomass(i,j)  * (1 - bc_fhydrophobic) * dt 
               qbrphobic_srcprm(i,kts,j) = f_smold(i,kts,j) *  mpas_emis_gocart2G%br_biomass(i,j)  *  br_fhydrophobic      * dt
               qbrphilic_srcprm(i,kts,j) = f_smold(i,kts,j) *  mpas_emis_gocart2G%br_biomass(i,j)  * (1 - br_fhydrophobic) * dt     
               qocphobic_srcprm(i,kts,j) = f_smold(i,kts,j) *  mpas_emis_gocart2G%oc_biomass(i,j)  *  oc_fhydrophobic      * dt
               qocphilic_srcprm(i,kts,j) = f_smold(i,kts,j) *  mpas_emis_gocart2G%oc_biomass(i,j)  * (1 - oc_fhydrophobic) * dt
               qso2_srcprm(i,kts,j)      = f_smold(i,kts,j) *  mpas_emis_gocart2G%su_biomass(i,j)     * dt 
               qnh3_srcprm(i,kts,j)      = f_smold(i,kts,j) *    mpas_emis_gocart2G%nh3_bb(i,j)       * dt
               qsoapbb_srcprm(i,kts,j)   = f_smold(i,kts,j) *  mpas_emis_gocart2G%soap_biomass(i,j)   * dt
               

               !Mixing ratio at PRM injH levels (flaming phase):
                                    
               do k = mpas_prm%k1_zbot(i,j), mpas_prm%k2_ztop(i,j)

                  f_flame(i,k,j) = mpas_prm%flame_fraction(i,j) / (mpas_prm%rho_dry(i,k,j) * mpas_prm%dz(i,k,j)) !multiplicar pelo dt aqui
                  
                  !The sum corrects the mixing ratio for the case when injH = 1       

                  qbcphobic_srcprm(i,k,j) = qbcphobic_srcprm(i,k,j) + (f_flame(i,k,j)  *  mpas_emis_gocart2G%bc_biomass(i,j)  *  bc_fhydrophobic      *  dt )
                  qbcphilic_srcprm(i,k,j) = qbcphilic_srcprm(i,k,j) + (f_flame(i,k,j)  *  mpas_emis_gocart2G%bc_biomass(i,j)  * (1 - bc_fhydrophobic) *  dt )
                  qbrphobic_srcprm(i,k,j) = qbrphobic_srcprm(i,k,j) + (f_flame(i,k,j)  *  mpas_emis_gocart2G%br_biomass(i,j)  *  br_fhydrophobic      *  dt )
                  qbrphilic_srcprm(i,k,j) = qbrphilic_srcprm(i,k,j) + (f_flame(i,k,j)  *  mpas_emis_gocart2G%br_biomass(i,j)  * (1 - br_fhydrophobic) *  dt ) 
                  qocphobic_srcprm(i,k,j) = qocphobic_srcprm(i,k,j) + (f_flame(i,k,j)  *  mpas_emis_gocart2G%oc_biomass(i,j)  *  oc_fhydrophobic      *  dt )
                  qocphilic_srcprm(i,k,j) = qocphilic_srcprm(i,k,j) + (f_flame(i,k,j)  *  mpas_emis_gocart2G%oc_biomass(i,j)  * (1 - oc_fhydrophobic) *  dt )
                  qso2_srcprm(i,k,j)      = qso2_srcprm(i,k,j)    + (f_flame(i,k,j)  *  mpas_emis_gocart2G%su_biomass(i,j)   *  dt )
                  qnh3_srcprm(i,k,j)      = qnh3_srcprm(i,k,j)    + (f_flame(i,k,j)  *     mpas_emis_gocart2G%nh3_bb(i,j)    *  dt )
                  qsoapbb_srcprm(i,k,j)   = qsoapbb_srcprm(i,k,j) + (f_flame(i,k,j)  *  mpas_emis_gocart2G%soap_biomass(i,j) *  dt )
                                                    
               enddo

            enddo
         enddo
      
      !----------------------------------------------------------------------------------
      ! Updates Wildfire Aerosol Mixing Ratio 
      !----------------------------------------------------------------------------------
         
         call mpas_pool_get_array(state,'scalars',scalars,time_lev)
         call mpas_pool_get_dimension(state,'index_qbcphobic',index_qbcphobic)
         call mpas_pool_get_dimension(state,'index_qbcphilic',index_qbcphilic)
         call mpas_pool_get_dimension(state,'index_qbrphobic',index_qbrphobic)
         call mpas_pool_get_dimension(state,'index_qbrphilic',index_qbrphilic)
         call mpas_pool_get_dimension(state,'index_qocphobic',index_qocphobic)
         call mpas_pool_get_dimension(state,'index_qocphilic',index_qocphilic)
         call mpas_pool_get_dimension(state,'index_qso2',index_qso2)
         call mpas_pool_get_dimension(state,'index_qnh3',index_qnh3)
         call mpas_pool_get_dimension(state,'index_qsoapbb',index_qsoapbb)

         do j = jts, jte
            do i = its, ite
               do k = kts, kte
                  
                  scalars(index_qbcphobic,k,i) = scalars(index_qbcphobic,k,i) + qbcphobic_srcprm(i,k,j) 
                  scalars(index_qbcphilic,k,i) = scalars(index_qbcphilic,k,i) + qbcphilic_srcprm(i,k,j)
                  scalars(index_qbrphobic,k,i) = scalars(index_qbrphobic,k,i) + qbrphobic_srcprm(i,k,j)
                  scalars(index_qbrphilic,k,i) = scalars(index_qbrphilic,k,i) + qbrphilic_srcprm(i,k,j)
                  scalars(index_qocphobic,k,i) = scalars(index_qocphobic,k,i) + qocphobic_srcprm(i,k,j)
                  scalars(index_qocphilic,k,i) = scalars(index_qocphilic,k,i) + qocphilic_srcprm(i,k,j)
                  scalars(index_qso2,k,i)    = scalars(index_qso2,k,i)     + qso2_srcprm(i,k,j)
                  scalars(index_qnh3,k,i)    = scalars(index_qnh3,k,i)     + qnh3_srcprm(i,k,j)
                  scalars(index_qsoapbb,k,i) = scalars(index_qsoapbb,k,i)  + qsoapbb_srcprm(i,k,j)

               enddo
            enddo
         enddo

         

      !-------------------------------------------------------------------------------------------------
      !Deallocates local arrays that received Mixing Ratio Tendency due to emissions
      !-------------------------------------------------------------------------------------------------
         
         if(allocated(f_smold))     deallocate(f_smold)
         if(allocated(f_flame))     deallocate(f_flame)

         if(allocated(qbcphobic_srcprm)) deallocate(qbcphobic_srcprm)
         if(allocated(qbcphilic_srcprm)) deallocate(qbcphilic_srcprm)
         if(allocated(qbrphobic_srcprm)) deallocate(qbrphobic_srcprm)
         if(allocated(qbrphilic_srcprm)) deallocate(qbrphilic_srcprm)
         if(allocated(qocphobic_srcprm)) deallocate(qocphobic_srcprm)
         if(allocated(qocphilic_srcprm)) deallocate(qocphilic_srcprm)
         if(allocated(qso2_srcprm))      deallocate(qso2_srcprm)
         if(allocated(qnh3_srcprm))      deallocate(qnh3_srcprm)
         if(allocated(qsoapbb_srcprm))   deallocate(qsoapbb_srcprm)


 call mpas_log_write('--- end subroutine gocart2G_plumerisemodel.')
 end subroutine gocart2G_plumerisemodel
!==================================================================================================================


!==================================================================================================================
 end module mpas_chemistry_gocart2G_plumerisemodel
!==================================================================================================================
