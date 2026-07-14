! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
! Deallocate Plume Rise Model Arrays 
! Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-06-14.
!
!==================================================================================================================
 module mpas_chemistry_plumerisemodel_finalize
 use mpas_log
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_config
 use mpas_chemistry_plumerisemodel_vars,only: mpas_prm, mpas_lbc_prm
 

 implicit none
 private
 public:: plumerisemodel_finalize


 contains


!==================================================================================================================
 subroutine plumerisemodel_finalize(configs)
!==================================================================================================================

    !*******************************************************************************
    !implicit none
    !*******************************************************************************

        !--- input arguments:
        type(mpas_pool_type),intent(in):: configs

        call mpas_log_write ('--- enter subroutine plumerisemodel_finalize:')

    !*******************************************************************************
    !Deallocate PRM Arrays
    !*******************************************************************************

        
        !--- Deallocate PRM Arrays (subroutine at mpas_chemistry_plumerisemodel_interface.F90)
        call mpas_prm%prm_deallocate()

        !--- Deallocate PRM lbc Arrays (subroutine at mpas_chemistry_plumerisemodel_lbc.F90)
        call mpas_lbc_prm%prm_deallocate()

 call mpas_log_write ('--- end subroutine plumerisemodel_finalize.')
 end subroutine plumerisemodel_finalize
!==================================================================================================================

!==================================================================================================================
end module mpas_chemistry_plumerisemodel_finalize
!==================================================================================================================

