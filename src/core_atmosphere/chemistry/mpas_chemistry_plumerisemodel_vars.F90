! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
!Defines derived-type objects to store Plume Rise Model related variables.
!Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-04-08.
!
!==================================================================================================================
 module mpas_chemistry_plumerisemodel_vars
 use mpas_chemistry_plumerisemodel_interface,only: atm_prm
 use mpas_chemistry_plumerisemodel_lbc,only: lbc_prm
 

    !*******************************************************************************
    !Declares PRM types
    !*******************************************************************************

        type(atm_prm):: mpas_prm
        type(lbc_prm):: mpas_lbc_prm
       

!==================================================================================================================
 end module mpas_chemistry_plumerisemodel_vars
!==================================================================================================================

