! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
!!Initialization of Plume Rise Model Input Variables (Fire Size and FRP).
!Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-04-08.
!
!==================================================================================================================
 module mpas_chemistry_plumerisemodel_lbc
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_array


 implicit none
 private

 type,public:: lbc_prm
    
    integer:: its,ite,jts,jte,kts,kte

    !--- biomass burning input arrays for PRM:
    real(kind=RKIND),dimension(:,:),pointer:: firesize_avg_biomass      => null()
    real(kind=RKIND),dimension(:,:),pointer:: firesize_std_biomass      => null()
    real(kind=RKIND),dimension(:,:),pointer:: frp_avg_biomass           => null()
    real(kind=RKIND),dimension(:,:),pointer:: frp_std_biomass           => null()

    contains
       procedure:: prm_allocate   => mpas_chemistry_plumerisemodel_lbc_allocate
       procedure:: prm_initlbc    => mpas_chemistry_plumerisemodel_initlbc

 end type

contains

!==================================================================================================================
 subroutine mpas_chemistry_plumerisemodel_lbc_allocate(self,its,ite,jts,jte,kts,kte)
!==================================================================================================================


    !*******************************************************************************
    !implicit none
    !*******************************************************************************

        !input arguments:
            integer,intent(in):: its,ite,jts,jte,kts,kte

        !inout arguments:
            class(lbc_prm),intent(inout):: self


        call mpas_log_write(' --- enter subroutine mpas_chemistry_plumerisemodel_lbc_allocate:')


    !*******************************************************************************
    !Allocate biomass burning input arrays for PRM:
    !*******************************************************************************

        
        if(.not.associated(self%firesize_avg_biomass) ) allocate(self%firesize_avg_biomass(its:ite,jts:jte)       )
        if(.not.associated(self%firesize_std_biomass) ) allocate(self%firesize_std_biomass(its:ite,jts:jte)       )
        if(.not.associated(self%frp_avg_biomass)     )  allocate(self%frp_avg_biomass(its:ite,jts:jte)            )
        if(.not.associated(self%frp_std_biomass)     )  allocate(self%frp_std_biomass(its:ite,jts:jte)            )



 call mpas_log_write(' --- end subroutine mpas_chemistry_plumerisemodel_lbc_allocate.')

 end subroutine mpas_chemistry_plumerisemodel_lbc_allocate
!==================================================================================================================


!==================================================================================================================
 subroutine mpas_chemistry_plumerisemodel_initlbc(self,biob_prm_lbc,its,ite,jts,jte,kts,kte)
!==================================================================================================================

    !*******************************************************************************
    !implicit none
    !*******************************************************************************
    

        !--- input arguments:
            type(mpas_pool_type),intent(in):: biob_prm_lbc
            integer,intent(in):: its,ite,jts,jte,kts,kte

        !--- inout arguments:
            class(lbc_prm),intent(inout):: self

        !--- local variables and pointers:  
            real(kind=RKIND),dimension(:),pointer:: firesize_avg, firesize_std, frp_avg, frp_std

        !--- local variables and arrays:
            integer:: i,j,k

        call mpas_log_write(' --- enter subroutine mpas_chemistry_plumerisemodel_lbcinit:')



    !*******************************************************************************
    !Fill local array with PRM input data initialization 
    !*******************************************************************************
    
        call mpas_pool_get_array(biob_prm_lbc,'firesize_avg',firesize_avg)
        call mpas_pool_get_array(biob_prm_lbc,'firesize_std',firesize_std)
        call mpas_pool_get_array(biob_prm_lbc,'frp_avg'     ,frp_avg     )
        call mpas_pool_get_array(biob_prm_lbc,'frp_std'     ,frp_std     )

        do j = jts,jte
            do i = its,ite
                self%firesize_avg_biomass(i,j) = firesize_avg(i)
                self%firesize_std_biomass(i,j) = firesize_std(i)
                self%frp_avg_biomass(i,j)      = frp_avg(i)
                self%frp_std_biomass(i,j)      = frp_std(i)
                
                !DEBUG - Tracking variables:
                    !call mpas_log_write('  ncell(i)=$i', intArgs=(/i/))
                    !call mpas_log_write('  firesize_avg(i)=$r', realArgs=(/firesize_avg(i)/))
                    !call mpas_log_write('  self%firesize_avg_biomass(i,j)=$r', realArgs=(/self%firesize_avg_biomass(i,j)/))
                    
            enddo
        enddo

        call mpas_log_write('     fill local arrays, ex:')
        call mpas_log_write('     ncell(i)=$i', intArgs=(/53/))
        call mpas_log_write('     firesize_avg(i)=$r', realArgs=(/firesize_avg(53)/))
        call mpas_log_write('     self%firesize_avg_biomass(i,j)=$r', realArgs=(/self%firesize_avg_biomass(53,1)/))
                    

 call mpas_log_write('--- end subroutine mpas_chemistry_plumerisemodel_initlbc.')

 end subroutine mpas_chemistry_plumerisemodel_initlbc
!==================================================================================================================





!==================================================================================================================
 end module mpas_chemistry_plumerisemodel_lbc
!==================================================================================================================

