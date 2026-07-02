!Initialization and Updates of Plume Rise Model Input Variables (Fire Size and FRP).
!Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-03-28.

!================================================================================================================== 
 module mpas_chemistry_plumerisemodel_init
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types,only: mpas_pool_type
 use mpas_chemistry_plumerisemodel_interface
 use mpas_chemistry_plumerisemodel_biob_input
 use mpas_chemistry_plumerisemodel_vars,only: mpas_prm, mpas_biob_input_prm

 
 contains

!==================================================================================================================
 subroutine plumerisemodel_init(mesh, prm_fire_input)
!==================================================================================================================


   !*******************************************************************************
    implicit none
   !*******************************************************************************

      !--- input arguments:
         type(mpas_pool_type),intent(in):: mesh
         type(mpas_pool_type),intent(in):: prm_fire_input
      
      !--- local variables:
         integer:: its,ite,jts,jte,kts,kte


      call mpas_log_write(' ')
      call mpas_log_write('--- enter subroutine plumerisemodel_init:')


   !*******************************************************************************
   !Initialization of dimensions and arrays needed to initialize PRM Input Data 
   !*******************************************************************************
   
      !--- initializes dimensions:
         
         call mpas_prm%prm_dims(mesh)
         
         its   = mpas_prm%its
         ite   = mpas_prm%ite
         jts   = mpas_prm%jts
         jte   = mpas_prm%jte
         kts   = mpas_prm%kts
         kte   = mpas_prm%kte
     

      !--- allocates local arrays needed to run PRM Input Data:
      call mpas_biob_input_prm%prm_allocate(its,ite,jts,jte,kts,kte)


      !--- initializes local arrays needed to run PRM Input Data emissions:
      call mpas_biob_input_prm%prm_emissions(prm_fire_input,its,ite,jts,jte,kts,kte)



 end subroutine plumerisemodel_init
!=================================================================================================================


!==================================================================================================================
 end module mpas_chemistry_plumerisemodel_init
!==================================================================================================================