! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
!Interface for MPAS - Plume Rise Model.
!Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-04-08.
!
!==================================================================================================================
 module mpas_chemistry_plumerisemodel_interface
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_config,mpas_pool_get_dimension
 

 implicit none
 private


 type,public:: atm_prm
    
    integer:: its,ite,jts,jte,kts,kte,ktep1

    contains
       procedure:: prm_dims       => mpas_chemistry_plumerisemodel_dims
       
 end type


contains

!==================================================================================================================
 subroutine mpas_chemistry_plumerisemodel_dims(self,mesh)
!==================================================================================================================

    !*******************************************************************************
    !implicit none
    !*******************************************************************************
        
        !--- input arguments:
            type(mpas_pool_type),intent(in):: mesh

        !--- inout arguments:
            class(atm_prm),intent(inout):: self

        !--- local variables:
            integer,pointer:: nCellsSolve,nVertLevels


        call mpas_log_write(' --- enter subroutine mpas_chemistry_plumerisemodel_dims:')

    !*******************************************************************************
    !Get MPAS mesh dimensions
    !*******************************************************************************
     
        call mpas_pool_get_dimension(mesh,'nCellsSolve'    ,nCellsSolve )
        call mpas_pool_get_dimension(mesh,'nVertLevels'    ,nVertLevels )

        self%its = 1 ; self%ite = nCellsSolve
        self%jts = 1 ; self%jte = 1
        self%kts = 1 ; self%kte = nVertLevels ; self%ktep1 = nVertLevels+1
        call mpas_log_write('     ITS = $i   ITE = $i',intArgs=(/self%its,self%ite/))
        call mpas_log_write('     JTS = $i   JTE = $i',intArgs=(/self%jts,self%jte/))
        call mpas_log_write('     KTS = $i   KTE = $i',intArgs=(/self%kts,self%kte/))

 
 call mpas_log_write(' --- end subroutine mpas_chemistry_plumerisemodel_dims.')

 end subroutine mpas_chemistry_plumerisemodel_dims
!==================================================================================================================




!==================================================================================================================
 end module mpas_chemistry_plumerisemodel_interface
!==================================================================================================================