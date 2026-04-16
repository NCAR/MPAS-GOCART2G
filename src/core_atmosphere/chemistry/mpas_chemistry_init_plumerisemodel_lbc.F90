! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
!Initialization of Plume Rise Model Input Variables (Fire Size and FRP).
!Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-03-31.
!
!================================================================================================================== 
module mpas_chemistry_init_plumerisemodel_lbc
use mpas_log
use mpas_kind_types
use mpas_pool_routines
use mpas_derived_types
use mpas_stream_manager

implicit none
 private
 public:: init_plumerisemodel_lbc
 
 contains


!==================================================================================================================
 subroutine init_plumerisemodel_lbc(clock,stream_manager,mesh,FINN_biob_prm_lbc, biob_prm_lbc)
!==================================================================================================================

   !*******************************************************************************
   !implicit none
   !*******************************************************************************

      !input arguments:
         type(mpas_Clock_type),intent(in),pointer:: clock
         type(mpas_pool_type), intent(in),pointer:: mesh
         type(mpas_pool_type), intent(in),pointer:: FINN_biob_prm_lbc
         

      !inout arguments:
         type(MPAS_streamManager_type),intent(inout):: stream_manager
         type(mpas_pool_type),intent(inout),pointer:: biob_prm_lbc


      !FINN Variables (Fire Size and FRP from pool FINN_biob_prm_lbc):
         real(kind=RKIND),dimension(:),pointer:: area_biob_modis
         real(kind=RKIND),dimension(:),pointer:: area_std_biob_modis
         real(kind=RKIND),dimension(:),pointer:: frp_biob_modis
         real(kind=RKIND),dimension(:),pointer:: frp_std_biob_modis

      !PRM Variables (Fire Size and FRP from pool biob_prm_lbc):
         real(kind=RKIND),dimension(:),pointer:: firesize_avg
         real(kind=RKIND),dimension(:),pointer:: firesize_std
         real(kind=RKIND),dimension(:),pointer:: frp_avg
         real(kind=RKIND),dimension(:),pointer:: frp_std

      !local variables and arrays:
         character(len=StrKIND):: actualTimeStamp
         integer:: iCell
         integer,pointer:: nCells
         

      call mpas_log_write(' ')
      call mpas_log_write('--- enter subroutine init_plumerisemodel_lbc:')


   !*******************************************************************************
   ! Initialization of PRM input variables - Fire Size and FRP
   !*******************************************************************************

      call mpas_pool_get_dimension(mesh,'nCells',nCells)
      call mpas_log_write('    nCells = $i',intArgs=(/nCells/))


      !----------------------------------------------------------------------
      !Initialize the PRM variables as zero
      !----------------------------------------------------------------------
         call mpas_pool_get_array(biob_prm_lbc,'firesize_avg', firesize_avg)
         call mpas_pool_get_array(biob_prm_lbc,'firesize_std', firesize_std)
         call mpas_pool_get_array(biob_prm_lbc,'frp_avg', frp_avg)
         call mpas_pool_get_array(biob_prm_lbc,'frp_std', frp_std)

            
         firesize_avg(:)  = 0._RKIND
         firesize_std(:)  = 0._RKIND
         frp_avg(:)       = 0._RKIND
         frp_std(:)       = 0._RKIND

      !----------------------------------------------------------------------
      !Read the FINN variables from input streams
      !----------------------------------------------------------------------
         
         
         !Read Active Fire Size Average
            call mpas_log_write('    read input stream: prm_lbc_area_avg')
            call mpas_pool_get_array(FINN_biob_prm_lbc,'area_biob_modis', area_biob_modis)
            call mpas_stream_mgr_read(stream_manager,'prm_lbc_area_avg',rightNow=.true., &
                     whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
            call mpas_log_write('    latest time before is  '//trim(actualTimestamp))
         
         !Read Active Fire Size Standard Deviation
            call mpas_log_write('    read input stream: prm_lbc_area_std')
            call mpas_pool_get_array(FINN_biob_prm_lbc,'area_std_biob_modis', area_std_biob_modis)
            call mpas_stream_mgr_read(stream_manager,'prm_lbc_area_std',rightNow=.true., &
                     whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
            
         !Read Fire Radiative Power Average
            call mpas_log_write('    read input stream: prm_lbc_frp_avg')
            call mpas_pool_get_array(FINN_biob_prm_lbc,'frp_biob_modis', frp_biob_modis)
            call mpas_stream_mgr_read(stream_manager,'prm_lbc_frp_avg',rightNow=.true., &
                     whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
            
         !Read Fire Radiative Power Standard Deviation
            call mpas_log_write('    read input stream: prm_lbc_frp_std')
            call mpas_pool_get_array(FINN_biob_prm_lbc,'frp_std_biob_modis', frp_std_biob_modis)
            call mpas_stream_mgr_read(stream_manager,'prm_lbc_frp_std',rightNow=.true., &
                     whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
            
         

         !DEBUG - Tracking variables:
         !do iCell = 1, nCells
         !   call mpas_log_write('  area_biob_modis(iCell)=$r', realArgs=(/area_biob_modis(iCell)/),messageType=MPAS_LOG_ERR)
         !enddo
         


         !--------------------------------------------------------------------
         !Initialize the PRM variables as FINN values
         !--------------------------------------------------------------------
                   
            firesize_avg(:) = area_biob_modis(:) 
            firesize_std(:) = area_std_biob_modis(:) 
            frp_avg(:)      = frp_biob_modis(:) 
            frp_std(:)      = frp_std_biob_modis(:) 



 call mpas_log_write('--- end subroutine init_plumerisemodel_lbc:')
 call mpas_log_write(' ')
 end subroutine init_plumerisemodel_lbc
!=================================================================================================================




!==================================================================================================================
 end module mpas_chemistry_init_plumerisemodel_lbc
!==================================================================================================================