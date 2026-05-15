! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module GA_EnvironmentMod
 use mpas_kind_types,only: RKIND
 use mpas_log

 use GOCART2G_instance,only: wavelengths_for_profile_aop_in_nm, &
                             wavelengths_for_vertically_integrated_aop_in_nm
 use GOCART2G_MieMod_smiol

 implicit none
 private
 public :: GA_Environment


 type:: GA_Environment
    type(GOCART2G_Mie):: rad_Mie,diag_Mie

!   logical:: scav_byColdCloud                                      ! new flag example

    integer:: rhFlag
    integer:: nbins
    integer:: km                                                    ! vertical grid dimension
    integer:: instance                                              ! data or computational instance
    integer:: klid                                                  ! vertical index of pressure lid

    real(kind=RKIND):: cdt                                          ! chemistry timestep (secs)
    real(kind=RKIND):: plid                                         ! pressure lid [hPa]
    real(kind=RKIND),dimension(:),allocatable:: radius              ! particle effective radius [um]
    real(kind=RKIND),dimension(:),allocatable:: rhop                ! soil class density [kg m-3]
    real(kind=RKIND),dimension(:),allocatable:: fscav               ! scavenging efficiency
    real(kind=RKIND),dimension(:),allocatable:: molwght             ! molecular weight
    real(kind=RKIND),dimension(:),allocatable:: fnum                ! number of particles per kg mass
    real(kind=RKIND),dimension(:),allocatable:: wavelengths_profile ! wavelengths for profile aop [nm]
    real(kind=RKIND),dimension(:),allocatable:: wavelengths_vertint ! wavelengths for vertically integrated aop [nm]

! Gas scavenging metadata (6 params per bin) ---
    real(kind=RKIND),dimension(:,:),allocatable:: scav_data_gas

    contains
       procedure:: load_from_config
 end type GA_Environment


!--- locals:
 integer:: status
 integer:: nbins
 integer:: n_wavelengths_profile,n_wavelengths_vertint,n_channels


 contains


!=================================================================================================================
 subroutine load_from_config(self,nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                             rhFlag,pressure_lid_in_hPa,scav_data_gas_in)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: nbins
 integer,intent(in):: rhFlag

 real(kind=RKIND),intent(in):: pressure_lid_in_hPa
 real(kind=RKIND),intent(in),dimension(nbins):: particle_radius_microns
 real(kind=RKIND),intent(in),dimension(nbins):: particle_density
 real(kind=RKIND),intent(in),dimension(nbins):: molecular_weight
 real(kind=RKIND),intent(in),dimension(nbins):: fscav
 real(kind=RKIND),intent(in),dimension(nbins):: fnum

!--- NEW: Optional input for Gas Scavenging ---
    real(kind=RKIND), intent(in), dimension(6, nbins), optional :: scav_data_gas_in

!--- inout arguments:
 class(GA_Environment),intent(inout) :: self


!--- local variables:
 integer:: n
 integer:: n_wavelengths_profile
 integer:: n_wavelengths_vertint

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
!call mpas_log_write('--- enter subroutine load_from_config:')

 self%nbins = nbins

!--- allocate arrays contained in GA_Environment:
 if(.not.allocated(self%radius) ) allocate(self%radius(nbins) )
 if(.not.allocated(self%rhop)   ) allocate(self%rhop(nbins)   )
 if(.not.allocated(self%fscav)  ) allocate(self%fscav(nbins)  )
 if(.not.allocated(self%molwght)) allocate(self%molwght(nbins))
 if(.not.allocated(self%fnum)   ) allocate(self%fnum(nbins)   )


!--- initialize arrays in GA_Environment:
 do n = 1, nbins
    self%radius(n)  = particle_radius_microns(n)
    self%rhop(n)    = particle_density(n)
    self%fscav(n)   = fscav(n)
    self%molwght(n) = molecular_weight(n)
    self%fnum(n)    = fnum(n)
 enddo

 self%rhFlag = rhFlag
 self%plid   = pressure_lid_in_hPa

!--- NEW: Initialize Gas Scavenging Data ---
    if(.not.allocated(self%scav_data_gas)) allocate(self%scav_data_gas(6, nbins))

    if (present(scav_data_gas_in)) then
        ! For NH3, SO2, etc., load the MOZCART constants
        self%scav_data_gas = scav_data_gas_in
    else
        ! For Dust, Sea Salt, etc., initialize to zero
        self%scav_data_gas = 0.0_RKIND
    end if

!--- DO NOT KNOW YET HOW TO INITIALIZE THOSE ARRAYS:
 n_wavelengths_profile = size(wavelengths_for_profile_aop_in_nm)
 n_wavelengths_vertint = size(wavelengths_for_vertically_integrated_aop_in_nm)
!call mpas_log_write('--- n_wavelengths_profile = $i',intArgs=(/n_wavelengths_profile/))
!call mpas_log_write('--- n_wavelengths_vertint = $i',intArgs=(/n_wavelengths_vertint/))

 if(.not.allocated(self%wavelengths_profile)) allocate(self%wavelengths_profile(n_wavelengths_profile))
 if(.not.allocated(self%wavelengths_vertint)) allocate(self%wavelengths_vertint(n_wavelengths_vertint))

 do n = 1, n_wavelengths_profile
    self%wavelengths_profile(n) = wavelengths_for_profile_aop_in_nm(n)
!   call mpas_log_write('$r',realArgs=(/self%wavelengths_profile(n)/))
 enddo
 do n = 1, n_wavelengths_vertint
    self%wavelengths_vertint(n) = wavelengths_for_vertically_integrated_aop_in_nm(n)
!   call mpas_log_write('$r',realArgs=(/self%wavelengths_vertint(n)/))
 enddo

!call mpas_log_write('--- end subroutine load_from_config.')

 end subroutine load_from_config

!=================================================================================================================
 end module GA_EnvironmentMod
!=================================================================================================================
