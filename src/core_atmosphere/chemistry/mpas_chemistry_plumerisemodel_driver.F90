! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
!Plume Rise Model Main Driver
!Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-04-20.
!
!==================================================================================================================
 module mpas_chemistry_plumerisemodel_driver

 use mpas_log
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines,only: mpas_pool_get_config,mpas_pool_get_subpool

 use mpas_chemistry_plumerisemodel_interface
 use mpas_chemistry_plumerisemodel_vars,only: mpas_prm, mpas_lbc_prm
 use mpas_chemistry_plumerisemodel_run,only: plumerisemodel_run


 implicit none
    
    private
    public:: plumerisemodel_driver

 contains

!==================================================================================================================
 subroutine plumerisemodel_driver(domain, mesh)
!==================================================================================================================


    !*******************************************************************************
    !implicit none
    !*******************************************************************************

        !inout arguments:

            type(domain_type),intent(inout):: domain

        !local variables and pointers:

            type(mpas_pool_type),pointer:: mesh
            type(mpas_pool_type),pointer:: state
            type(mpas_pool_type),pointer:: diag
            type(mpas_pool_type),pointer:: sfc_input
            type(mpas_pool_type),pointer:: diag_prm
            type(block_type),pointer:: block

            integer:: time_lev
            integer:: i,its,ite,j,jts,jte,k,kts,kte

  
        call mpas_log_write('--- enter subroutine plumerisemodel_driver:')

    !*******************************************************************************
    ! PRM Setup and Run Procedure
    !*******************************************************************************
        
        time_lev = 2

        block => domain % blocklist
        do while(associated(block))

            call mpas_pool_get_subpool(block%structs,'mesh'                ,mesh                )
            call mpas_pool_get_subpool(block%structs,'state'               ,state               )
            call mpas_pool_get_subpool(block%structs,'diag'                ,diag                )
            call mpas_pool_get_subpool(block%structs,'sfc_input'           ,sfc_input           )
            call mpas_pool_get_subpool(block%structs,'diag_prm'            ,diag_prm            )

            !-------------------------------------------------------------------------------
            ! Defines dimensions needed to run PRM:
            !-------------------------------------------------------------------------------    
                call mpas_prm%prm_dims(mesh)
                
                its    = mpas_prm%its
                ite    = mpas_prm%ite
                jts    = mpas_prm%jts
                jte    = mpas_prm%jte
                kts    = mpas_prm%kts
                kte    = mpas_prm%kte

            !-------------------------------------------------------------------------------
            ! Allocate local arrays from MPAS needed to run PRM:
            !-------------------------------------------------------------------------------
                
                call mpas_prm%prm_allocate()

            !-------------------------------------------------------------------------------
            ! Fills local chemistry arrays with global chemistry arrays:
            !-------------------------------------------------------------------------------

                call mpas_prm%prm_fromMPAS(mesh, state, diag, sfc_input, time_lev)

           
            !-------------------------------------------------------------------------------
            ! Run PRM 
            !-------------------------------------------------------------------------------

                !--- Determines the flame fraction and  combined vegetation from MPAS to PRM categories 

                    !-----------------------------------------------------------------------------------------------
                    !Flaming Consumption Freitas et al. (2006)  | LANDUSE.TBL = MODIFIED IGBP MODIS NOAH Land Cover! 
                    !-----------------------------------------------------------------------------------------------
                    ! 45,   ! Tropical Forest        |  igbp 1 thru 5 (Evergreen|Deciduos Needleleaf|Broadleaf Forest and Mixed)
                    ! 75,   ! cerrado/woody savanna  |  igbp 6 thru 9 and 11 (Shrublands, savannas, Permanent wetlands(pantanal)?)
                    ! 97    ! Grassland/cropland     |  igbp 10, 12 and 14 (grassland, cropland, natural cropland)
                    !  0    ! case default           !  Tundra, barren vegetated, water, urban, ice
                    !----------------------------------------------------------------------------------------------


                    do j = jts,jte
                        do i = its, ite
                            select case(mpas_prm%ivgtyp(i, j))
                                
                                case(2)
                                    mpas_prm%iveg_ag(i,j) = 1                   !Tropical Forest        !igbp = Evergreen Broadleaf Forest
                                    mpas_prm%flame_fraction(i,j) = 0.45_RKIND

                                case(1,3:5)
                                    mpas_prm%iveg_ag(i,j) = 2                   !Boreal Forest          !igbp = Evergreen|Deciduos Needleleaf, Deciduous Broadleaf, and Mixed
                                    mpas_prm%flame_fraction(i,j) = 0.45_RKIND
                                
                                case(6:9,11)
                                    mpas_prm%iveg_ag(i,j) = 3                   !Cerrado/woody savanna  !igbp = Shrublands, savannas, Permanent wetlands(pantanal)?
                                    mpas_prm%flame_fraction(i,j) = 0.75_RKIND
                                
                                case(10,12,14)
                                    mpas_prm%iveg_ag(i,j) = 4                   !Grassland/cropland     !igbp = grassland, cropland, natural cropland
                                    mpas_prm%flame_fraction(i,j) = 0.97_RKIND
                                
                                case default
                                    mpas_prm%iveg_ag(i,j) = 0                   ! Default value if not in the range
                                    mpas_prm%flame_fraction(i,j) = 0._RKIND

                            end select

                        enddo
                    enddo
                    
                
                !--- Corrects the PRM data input (Active Fire Size, FRP)

                    !There is some inconsistence between the wildfire data (emissions, firesize, frp)
                    !and the land use classification. Because of the interpolation from a lat lon 
                    !gridded file to the MPAS grid, some coastal areas that have fires, also has 
                    !the dominant vegetation type as ocean (because of the grid cell over the land
                    !and ocean). Therefore, we need to set the emissions over the ocean as zero. 
                    !For MODIFIED IGBP MODIS NOAH Land Cover: (ocean: ivgtyp=17)
                    !Since this behavior can also happen for other land use interfaces, 
                    !we corrected the data using the iveg_ag =0. 


                    do j = jts,jte
                        do i = its, ite

                            if (mpas_prm%iveg_ag(i,j) .eq. 0) then
                                mpas_prm%firesize_avg_corr(i,j)  = 0._RKIND
                                mpas_prm%firesize_std_corr(i,j)  = 0._RKIND
                                mpas_prm%frp_avg_corr(i,j)       = 0._RKIND
                                mpas_prm%frp_std_corr(i,j)       = 0._RKIND
                            else
                                mpas_prm%firesize_avg_corr(i,j)  = mpas_lbc_prm%firesize_avg_biomass(i,j)
                                mpas_prm%firesize_std_corr(i,j)  = mpas_lbc_prm%firesize_std_biomass(i,j)
                                mpas_prm%frp_avg_corr(i,j)       = mpas_lbc_prm%frp_avg_biomass(i,j)
                                mpas_prm%frp_std_corr(i,j)       = mpas_lbc_prm%frp_std_biomass(i,j)
                            endif

                        enddo
                    enddo

                     
                !--- Run the PRM to calculate smoke injection height

                    call plumerisemodel_run( )

                    
            !-------------------------------------------------------------------------------
            !Fills global chemistry arrays with local chemistry arrays:
            !-------------------------------------------------------------------------------
                
                call mpas_prm%prm_toMPAS(diag_prm)


            !-------------------------------------------------------------------------------
            ! Deallocate local arrays: (need to move to after gocart2 run)
            !-------------------------------------------------------------------------------

                call mpas_prm%prm_deallocate()

        block => block % next
        end do

 
 call mpas_log_write('--- end subroutine plumerisemodel_driver.')
 end subroutine plumerisemodel_driver
!==================================================================================================================


!==================================================================================================================
 end module mpas_chemistry_plumerisemodel_driver
!==================================================================================================================
