! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
!Run Plume Rise Model to Calculate Smoke Plume Injection height
!Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-04-20.
!
!==================================================================================================================
 module mpas_chemistry_plumerisemodel_run

 use mpas_log
 use mpas_derived_types 
 use mpas_kind_types,only: kind_phys => RKIND

 implicit none
    
    private                      
    public:: plumerisemodel_run



 contains

!==================================================================================================================
 subroutine plumerisemodel_run
!==================================================================================================================

 call mpas_log_write(' --- enter subroutine plumerisemodel_run:')




 call mpas_log_write(' --- end subroutine plumerisemodel_run.')
 end subroutine plumerisemodel_run
!==================================================================================================================


!==================================================================================================================
 end module mpas_chemistry_plumerisemodel_run
!==================================================================================================================
