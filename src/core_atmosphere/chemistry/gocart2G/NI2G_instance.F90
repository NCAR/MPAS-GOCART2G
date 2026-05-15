! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module NI2G_instance
 use mpas_kind_types,only: RKIND
 public
 save

!this module is the resource file for nitrate parameters. it is the same as NI2G_instance_NI.rc in the GOCART-2G
!directory ./GOCART-2G/ESMF/GOCART2G_GridComp/NI2G_GridComp.

!==================================================================================================================

!the location of those two files needs to be added if used:
!aerosol_radBands_optics_file: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_NI.v2_5.RRTMG.nc
!aerosol_monochromatic_optics_file: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_NI.v2_5.nc

!number of bins:
 integer,parameter:: nbins = 5

!scavenging efficiency per bin [km-1]:
 real(kind=RKIND),dimension(nbins),parameter:: fscav = (/0.0,0.4,0.4,0.4,0.4/)

!dry particle radius [um], used for settling:
 real(kind=RKIND),dimension(nbins),parameter:: particle_radius_microns =  (/0.0,0.2695,0.2695,2.1,7.57/)

!type of settling to use (see Chem_SettlingMod):
 integer,parameter:: rhFlag = 0

!dry particle density [kg m-3]:
 real(kind=RKIND),dimension(nbins),parameter:: particle_density =(/1000,1769,1725,2200,2650/)

 real(kind=RKIND),parameter:: pressure_lid_in_hPa = 0.01

!molecular weight of species [kg mole-1]:
 real(kind=RKIND),dimension(nbins),parameter:: molecular_weight = (/0.18,0.18,0.18,0.18,0.18/)

!number of particles per kg mass:
 real(kind=RKIND),dimension(nbins),parameter:: fnum = (/1.50e19,1.50e19,1.50e19,1.50e19,1.50e19/)

!number median radius [um]:
 real(kind=RKIND),dimension(nbins),parameter:: particle_radius_number = (/0.0118,0.0118,0.0118,0.0118,0.0118/)

!sigma of lognormal number distribution:
 real(kind=RKIND),dimension(nbins),parameter:: sigma = (/2.0,2.0,2.0,2.0,2.0/)

! Add to the parameter section
! Parameters: [H298 (M/atm), dHR (K), K1 (Kb), fK1 (dfKb), K2 (Kw), fK2 (dfKw)]
 real(kind=RKIND), dimension(6, nbins), parameter :: scav_data_gas = reshape((/ &
    60.2, 4160.0, 1.7e-5, -4325.0, 1.0e-14, -6716.0,  & ! Bin 1: NH3 (Gas)
    0.0,  0.0,    0.0,     0.0,     0.0,     0.0,     & ! Bin 2: NH4a
    0.0,  0.0,    0.0,     0.0,     0.0,     0.0,     & ! Bin 3: NO3an1
    0.0,  0.0,    0.0,     0.0,     0.0,     0.0,     & ! Bin 4: NO3an2
    0.0,  0.0,    0.0,     0.0,     0.0,     0.0      & ! Bin 5: NO3an3
 /), (/6, nbins/))

!==================================================================================================================
 end module NI2G_instance
!==================================================================================================================
