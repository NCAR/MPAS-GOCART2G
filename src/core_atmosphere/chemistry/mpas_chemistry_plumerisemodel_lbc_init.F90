! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
!Initialization of Plume Rise Model Input Variables (Fire Size and FRP).
!Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-04-08.
!
!================================================================================================================== 
 module mpas_chemistry_plumerisemodel_lbc_init
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types,only: mpas_pool_type
 use mpas_chemistry_plumerisemodel_vars, only: mpas_prm, mpas_lbc_prm
 use mpas_chemistry_plumerisemodel_interface
 use mpas_chemistry_plumerisemodel_lbc
 
 contains

!==================================================================================================================
 subroutine plumerisemodel_lbc_init(mesh, biob_prm_lbc)
!==================================================================================================================

    !*******************************************************************************
    ! implicit none
    !*******************************************************************************
        
        !--- input arguments:
            type(mpas_pool_type),intent(in):: mesh
            type(mpas_pool_type),intent(in):: biob_prm_lbc

        !--- local variables:
            integer:: its,ite,jts,jte,kts,kte


        call mpas_log_write(' ')
        call mpas_log_write('--- enter subroutine plumerisemodel_lbc_init:')


    !*******************************************************************************
    ! Initializes dimensions needed to initialize PRM input data:
    !*******************************************************************************

        call mpas_prm%prm_dims(mesh)

        its   = mpas_prm%its
        ite   = mpas_prm%ite
        jts   = mpas_prm%jts
        jte   = mpas_prm%jte
        kts   = mpas_prm%kts
        kte   = mpas_prm%kte

    !*******************************************************************************
    ! Allocates local arrays needed to run Plume Rise Model:
    !*******************************************************************************

        call mpas_lbc_prm%prm_allocate(its,ite,jts,jte,kts,kte)


    
    !*******************************************************************************
    ! Initializes local arrays needed to run Plume Rise Model:
    !*******************************************************************************
    
        call mpas_lbc_prm%prm_initlbc(biob_prm_lbc,its,ite,jts,jte,kts,kte)



 call mpas_log_write('--- end subroutine plumerisemodel_lbc_init.')
 call mpas_log_write(' ')
 end subroutine plumerisemodel_lbc_init
!=================================================================================================================




!==================================================================================================================
 end module mpas_chemistry_plumerisemodel_lbc_init
!==================================================================================================================