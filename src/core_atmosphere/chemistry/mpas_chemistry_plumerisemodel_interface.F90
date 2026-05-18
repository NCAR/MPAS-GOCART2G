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
 use mpas_atmphys_constants,only: R_d, R_v
  

 implicit none
 private


 type,public:: atm_prm
    
    integer:: its,ite,jts,jte,kts,kte,ktep1

    !--- atmospheric fields (from MPAS):
        
        real(kind=RKIND),dimension(:,:,:),pointer:: z                   => null()
        real(kind=RKIND),dimension(:,:,:),pointer:: zmid                => null()
        real(kind=RKIND),dimension(:,:,:),pointer:: dz                  => null()
        real(kind=RKIND),dimension(:,:,:),pointer:: t                   => null()
        real(kind=RKIND),dimension(:,:,:),pointer:: th                  => null()
        real(kind=RKIND),dimension(:,:,:),pointer:: u                   => null()
        real(kind=RKIND),dimension(:,:,:),pointer:: v                   => null()
        real(kind=RKIND),dimension(:,:,:),pointer:: rho_moist           => null()
        real(kind=RKIND),dimension(:,:,:),pointer:: pres                => null()
        real(kind=RKIND),dimension(:,:,:),pointer:: qvapor              => null()
        
        integer        , dimension(:,:)  ,pointer:: ivgtyp              => null()


    !--- smoke plume fields (from PRM):

        integer        , dimension(:,:)  ,pointer:: iveg_ag             => null()
        real(kind=RKIND),dimension(:,:)  ,pointer:: flame_fraction      => null()

        real(kind=RKIND),dimension(:,:)  ,pointer:: firesize_avg_corr   => null()
        real(kind=RKIND),dimension(:,:)  ,pointer:: firesize_std_corr   => null()
        real(kind=RKIND),dimension(:,:)  ,pointer:: frp_avg_corr        => null()
        real(kind=RKIND),dimension(:,:)  ,pointer:: frp_std_corr        => null()

        integer        , dimension(:,:)  ,pointer:: k1_zbot             => null()
        integer        , dimension(:,:)  ,pointer:: k2_ztop             => null()
        real(kind=RKIND),dimension(:,:)  ,pointer:: zbot                => null()
        real(kind=RKIND),dimension(:,:)  ,pointer:: ztop                => null()

        

    contains
        procedure:: prm_dims       => mpas_chemistry_plumerisemodel_dims
        procedure:: prm_allocate   => mpas_chemistry_plumerisemodel_allocate
        procedure:: prm_fromMPAS   => mpas_chemistry_plumerisemodel_fromMPAS
        procedure:: prm_toMPAS     => mpas_chemistry_plumerisemodel_toMPAS
        procedure:: prm_deallocate => mpas_chemistry_plumerisemodel_deallocate
       
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


        !call mpas_log_write(' --- enter subroutine mpas_chemistry_plumerisemodel_dims:')

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

 
 !call mpas_log_write(' --- end subroutine mpas_chemistry_plumerisemodel_dims.')

 end subroutine mpas_chemistry_plumerisemodel_dims
!==================================================================================================================



!==================================================================================================================
 subroutine mpas_chemistry_plumerisemodel_allocate(self)
!==================================================================================================================

    !*******************************************************************************
    !implicit none
    !*******************************************************************************

        !inout arguments:
            class(atm_prm),intent(inout):: self

        !local variables and arrays:
            integer:: its,ite,jts,jte,kts,kte,ktep1
  

      
        !call mpas_log_write(' --- enter subroutine mpas_chemistry_plumerisemodel_allocate:')

    !*******************************************************************************
    !Allocate local arrays needed to run PRM
    !*******************************************************************************


        its   = self%its
        ite   = self%ite
        jts   = self%jts
        jte   = self%jte
        kts   = self%kts
        kte   = self%kte

        
        !--- allocate atmospheric fields:

            if(.not.associated(self%z)    )             allocate( self%z                   (   its:ite, kts:kte, jts:jte)  )
            if(.not.associated(self%zmid) )             allocate( self%zmid                (   its:ite, kts:kte, jts:jte)  )
            if(.not.associated(self%dz)   )             allocate( self%dz                  (   its:ite, kts:kte, jts:jte)  )
            if(.not.associated(self%th)   )             allocate( self%th                  (   its:ite, kts:kte, jts:jte)  )
            if(.not.associated(self%t)    )             allocate( self%t                   (   its:ite, kts:kte, jts:jte)  )        
            if(.not.associated(self%u)    )             allocate( self%u                   (   its:ite, kts:kte, jts:jte)  )
            if(.not.associated(self%v)    )             allocate( self%v                   (   its:ite, kts:kte, jts:jte)  )
            if(.not.associated(self%rho_moist))         allocate( self%rho_moist           (   its:ite, kts:kte, jts:jte)  )
            if(.not.associated(self%pres))              allocate( self%pres                (   its:ite, kts:kte, jts:jte)  )
            if(.not.associated(self%qvapor))            allocate( self%qvapor              (   its:ite, kts:kte, jts:jte)  )
            if(.not.associated(self%ivgtyp))            allocate( self%ivgtyp              (   its:ite, jts:jte)           )


        !--- allocate smoke plume fields:

            if(.not.associated(self%iveg_ag))           allocate( self%iveg_ag             (   its:ite, jts:jte)           )
            if(.not.associated(self%flame_fraction))    allocate( self%flame_fraction      (   its:ite, jts:jte)           )

            if(.not.associated(self%firesize_avg_corr)) allocate( self%firesize_avg_corr   (   its:ite, jts:jte)           )
            if(.not.associated(self%firesize_std_corr)) allocate( self%firesize_std_corr   (   its:ite, jts:jte)           )
            if(.not.associated(self%frp_avg_corr))      allocate( self%frp_avg_corr        (   its:ite, jts:jte)           )
            if(.not.associated(self%frp_std_corr))      allocate( self%frp_std_corr        (   its:ite, jts:jte)           )

            if(.not.associated(self%k1_zbot))           allocate( self%k1_zbot             (   its:ite, jts:jte)           )
            if(.not.associated(self%k2_ztop))           allocate( self%k2_ztop             (   its:ite, jts:jte)           )
            if(.not.associated(self%zbot))              allocate( self%zbot                (   its:ite, jts:jte)           )
            if(.not.associated(self%ztop))              allocate( self%ztop                (   its:ite, jts:jte)           )

            

 !call mpas_log_write(' --- end subroutine mpas_chemistry_plumerisemodel_allocate.')
 end subroutine mpas_chemistry_plumerisemodel_allocate
!==================================================================================================================


!==================================================================================================================
 subroutine mpas_chemistry_plumerisemodel_fromMPAS(self, mesh, state, diag, sfc_input, time_lev)
!==================================================================================================================

    !*******************************************************************************
    !implicit none
    !*******************************************************************************

        !--- input arguments:
            type(mpas_pool_type),intent(in):: mesh
            type(mpas_pool_type),intent(in):: state
            type(mpas_pool_type),intent(in):: diag
            type(mpas_pool_type),intent(in):: sfc_input
            integer,intent(in):: time_lev

        !--- inout arguments:
            class(atm_prm),intent(inout):: self

        !--- local variables and pointers for atmospheric fields:
            real(kind=RKIND),dimension(:,:),pointer:: zgrid
            real(kind=RKIND),dimension(:,:),pointer:: theta_m
            real(kind=RKIND),dimension(:,:),pointer:: exner
            real(kind=RKIND),dimension(:,:),pointer:: u,v
            real(kind=RKIND),dimension(:,:),pointer:: rho
            real(kind=RKIND),dimension(:,:),pointer:: pressure_p, pressure_b
            real(kind=RKIND),dimension(:,:,:),pointer:: scalars
            integer,pointer:: index_qv
            real(kind=RKIND),dimension(:,:),pointer:: qv
           
            integer,dimension(:),pointer:: ivgtyp
            
            integer:: i,its,ite,j,jts,jte,k,kts,kte


        
        !call mpas_log_write(' --- enter subroutine mpas_chemistry_plumerisemodel_fromMPAS:')


    !*******************************************************************************
    !Fills local arrays (i,j,k) with global arrays (i,k)
    !*******************************************************************************

        !-------------------------------------------------------------------------------
        ! Get dimensions
        !-------------------------------------------------------------------------------
            
            its   = self%its
            ite   = self%ite
            jts   = self%jts
            jte   = self%jte
            kts   = self%kts
            kte   = self%kte

        !-------------------------------------------------------------------------------
        ! Fills atmospheric fields:
        !-------------------------------------------------------------------------------

            call mpas_pool_get_array(mesh,      'zgrid'                 ,  zgrid                )
            call mpas_pool_get_array(state,     'theta_m'               ,  theta_m  ,   time_lev)
            call mpas_pool_get_array(diag,      'exner'                 ,  exner                )
            call mpas_pool_get_array(diag,      'uReconstructZonal'     ,  u                    )
            call mpas_pool_get_array(diag,      'uReconstructMeridional',  v                    )
            call mpas_pool_get_array(diag,      'rho'                   ,  rho                  ) 
            call mpas_pool_get_array(diag,      'pressure_p'            ,  pressure_p           )
            call mpas_pool_get_array(diag,      'pressure_base'         ,  pressure_b           )
            call mpas_pool_get_array(sfc_input, 'ivgtyp'                ,  ivgtyp               )

            call mpas_pool_get_dimension(state  ,'index_qv'             ,  index_qv             )
            call mpas_pool_get_array(state      ,'scalars'              ,  scalars,     time_lev)
            qv => scalars(index_qv,:,:)
            

            
            do j = jts, jte
                do k = kts, kte
                    do i = its, ite

                        !--- arrays located at w points:
                        self%z(i,k,j)           = zgrid(k,i)

                        !--- arrays located at theta points:
                        self%qvapor(i,k,j)      = max(0._RKIND,qv(k,i))

                        self%zmid(i,k,j)        = 0.5 * (zgrid(k+1,i) + zgrid(k,i))
                        self%dz(i,k,j)          = zgrid(k+1,i) - zgrid(k,i) 

                        self%th(i,k,j)          = theta_m(k,i) / (1._RKIND + R_v/R_d * max(0._RKIND,qv(k,i))) 
                        self%t(i,k,j)           = self%th(i,k,j)*exner(k,i)

                        self%u(i,k,j)           = u(k,i) 
                        self%v(i,k,j)           = v(k,i) 

                        self%rho_moist(i,k,j)   = rho(k,i)*(1._RKIND + max(0._RKIND,qv(k,i)) )
                        self%pres(i,k,j)        = pressure_p(k,i) + pressure_b(k,i)

                    enddo
                enddo
            enddo

            
            do j = jts,jte
                do i = its,ite
                    self%ivgtyp(i,j)   = ivgtyp(i)
                enddo
            enddo


        !-------------------------------------------------------------------------------
        ! Initialize PRM fields 
        !-------------------------------------------------------------------------------

            !--- Initialize PRM input fields 

                do j = jts,jte
                    do i = its,ite

                        self%iveg_ag(i,j)               =  0
                        self%flame_fraction(i,j)        =  0._RKIND

                        self%firesize_avg_corr(i,j)     =  0._RKIND
                        self%firesize_std_corr(i,j)     =  0._RKIND
                        self%frp_avg_corr(i,j)          =  0._RKIND
                        self%frp_std_corr(i,j)          =  0._RKIND

                    enddo
                enddo


            !--- Initialize PRM output fields
                do j = jts,jte
                    do i = its,ite
                    
                        ! Setup to first level the local arrays to receive index level for lower and upper plume injection heights
                        self%k1_zbot(i,j)  = 1
                        self%k2_ztop(i,j)  = 1
                        
                        ! Setup to zero the local arrays to receive lower and upper plume injection heights
                        self%zbot(i,j)  = 0._RKIND
                        self%ztop(i,j)  = 0._RKIND

                    enddo
                enddo
            

 !call mpas_log_write(' --- end subroutine mpas_chemistry_plumerisemodel_fromMPAS.')
 end subroutine mpas_chemistry_plumerisemodel_fromMPAS
!==================================================================================================================

 
!==================================================================================================================
 subroutine mpas_chemistry_plumerisemodel_toMPAS(self, diag_prm)
!==================================================================================================================

    !*******************************************************************************
    !implicit none
    !*******************************************************************************
        
        !--- input arguments:
            type(mpas_pool_type),intent(in):: diag_prm

        !--- inout arguments:
            class(atm_prm),intent(inout):: self

        !--- local variables and pointers:

            integer,dimension(:),pointer:: k1_injH_bot, k2_injH_top
            real(kind=RKIND), dimension(:), pointer:: injH_bot, injH_top
            integer:: i,its,ite,j,jts,jte,k,kts,kte
       


        !call mpas_log_write(' --- enter subroutine mpas_chemistry_plumerisemodel_toMPAS:')


        !*******************************************************************************
        !Fills global PRM fields
        !*******************************************************************************
        

        !--- Get Dimensions
            its = self%its
            ite = self%ite
            jts = self%jts
            jte = self%jte
            kts = self%kts
            kte = self%kte

        !--- Fills global PRM arrays
           
            call mpas_pool_get_array(diag_prm,'k1_injH_bot', k1_injH_bot )
            call mpas_pool_get_array(diag_prm,'k2_injH_top', k2_injH_top )
            call mpas_pool_get_array(diag_prm,'injH_bot', injH_bot)
            call mpas_pool_get_array(diag_prm,'injH_top', injH_top )
            
            do j = jts,jte
                do i = its,ite
                    k1_injH_bot(i)  = self%k1_zbot(i,j)
                    k2_injH_top(i)  = self%k2_ztop(i,j)
                    injH_bot(i)     = self%zbot(i,j)
                    injH_top(i)     = self%ztop(i,j)
                enddo
            enddo
  


 !call mpas_log_write(' --- end subroutine mpas_chemistry_plumerisemodel_toMPAS.')
 end subroutine mpas_chemistry_plumerisemodel_toMPAS
!==================================================================================================================



!==================================================================================================================
 subroutine mpas_chemistry_plumerisemodel_deallocate(self)
!==================================================================================================================

    !*******************************************************************************
    !implicit none
    !*******************************************************************************

        !--- inout arguments:
            class(atm_prm),intent(inout):: self

        !call mpas_log_write(' --- enter subroutine mpas_chemistry_plumerisemodel_deallocate:')


    !*******************************************************************************
    !Deallocate local arrays 
    !*******************************************************************************

        
        !--- allocate atmospheric fields:

            if(associated(self%z)    )              deallocate( self%z                  )
            if(associated(self%zmid) )              deallocate( self%zmid               ) 
            if(associated(self%dz)   )              deallocate( self%dz                 )
            if(associated(self%th)   )              deallocate( self%th                 )
            if(associated(self%t)    )              deallocate( self%t                  )        
            if(associated(self%u)    )              deallocate( self%u                  )
            if(associated(self%v)    )              deallocate( self%v                  )
            if(associated(self%rho_moist))          deallocate( self%rho_moist          )
            if(associated(self%pres))               deallocate( self%pres               )
            if(associated(self%qvapor))             deallocate( self%qvapor             )
            if(associated(self%ivgtyp))             deallocate( self%ivgtyp             )


        !--- allocate smoke plume fields:

            if(associated(self%k1_zbot))            deallocate( self%k1_zbot            )
            if(associated(self%k2_ztop))            deallocate( self%k2_ztop            )
            if(associated(self%zbot))               deallocate( self%zbot               )
            if(associated(self%ztop))               deallocate( self%ztop               )  

            if(associated(self%iveg_ag))            deallocate( self%iveg_ag            )  
            if(associated(self%flame_fraction))     deallocate( self%flame_fraction     ) 

            if(associated(self%firesize_avg_corr))  deallocate( self%firesize_avg_corr  )  
            if(associated(self%firesize_std_corr))  deallocate( self%firesize_std_corr  ) 
            if(associated(self%frp_avg_corr))       deallocate( self%frp_avg_corr       ) 
            if(associated(self%frp_std_corr))       deallocate( self%frp_std_corr       ) 


  !call mpas_log_write(' --- end subroutine mpas_chemistry_plumerisemodel_deallocate.')
 end subroutine mpas_chemistry_plumerisemodel_deallocate
!==================================================================================================================

!==================================================================================================================
 end module mpas_chemistry_plumerisemodel_interface
!==================================================================================================================