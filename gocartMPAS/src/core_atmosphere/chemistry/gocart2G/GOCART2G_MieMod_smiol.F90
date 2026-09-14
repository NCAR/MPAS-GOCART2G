! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module gocart2G_MieMod_smiol
 use mpas_kind_types, only: R4KIND,RKIND
 use mpas_derived_types
 use mpas_log,only: mpas_log_write

 implicit none
 private
 public:: GOCART2G_Mie

 integer,parameter:: NRH_BINS = 991

 type GOCART2G_Mie
     character(len=:),allocatable :: table_name
     integer:: nch  ! number of channels in table (replacement of nlamfda)
     integer:: nrh  ! number of RH values in table
     integer:: nbin ! number of size bins in table
     integer:: nMom ! number of moments of phase function
     integer:: nPol ! number of elements of scattering phase matrix

     !--- pointers available in all netCDF files:
     real(kind=RKIND),dimension(:),pointer        :: rh      => null() ! (r) RH values   [fraction]
     real(kind=RKIND),dimension(:,:),pointer      :: reff    => null() ! (r,b) effective radius [m]
     real(kind=RKIND),dimension(:,:,:),pointer    :: bext    => null() ! (r,c,b) bext values [m2 kg-1]
     real(kind=RKIND),dimension(:,:,:),pointer    :: bsca    => null() ! (r,c,b) bsca values [m2 kg-1]
     real(kind=RKIND),dimension(:,:,:),pointer    :: bbck    => null() ! (r,c,b) bbck values [m2 kg-1]
     real(kind=RKIND),dimension(:,:,:),pointer    :: g       => null() ! (r,c,b) asymmetry parameter
     real(kind=RKIND),dimension(:,:,:),pointer    :: refr    => null() ! (r,c,b) real part of refractive index
     real(kind=RKIND),dimension(:,:,:),pointer    :: refi    => null() ! (r,c,b) imaginary part of refractive index

     !--- pointers available in netCDF files with the "wavelength" option:
     real(kind=RKIND),dimension(:),pointer        :: wavelengths => null() ! (c) wavelengths [m]
     real(kind=RKIND),dimension(:,:,:,:),pointer  :: pback  => null() ! (r,c,b,m,p) backscatter phase function
     
     !20240601_rkumar_start
     ! real(kind=RKIND),dimension(:,:,:,:,:),pointer:: pmom    => null() ! (r,c,b,m,p) moments of phase function
     !20240601_rkumar_end

     !--- pointers (and derived pointers) that are sometimes available in netCDF with the "wavelength" option:
     real(kind=RKIND),dimension(:,:),pointer      :: gf      => null() ! (r,b) hygroscopic growth factor
     real(kind=RKIND),dimension(:,:),pointer      :: rhop    => null() ! (r,b) wet particle density [kg m-3]
     real(kind=RKIND),dimension(:,:),pointer      :: rhod    => null() ! (r,b) wet particle density [kg m-3]
     real(kind=RKIND),dimension(:,:),pointer      :: vol     => null() ! (r,b) wet particle volume [m3 kg-1]
     real(kind=RKIND),dimension(:,:),pointer      :: area    => null() ! (r,b) wet particle cross section [m2 kg-1]

     real(kind=RKIND),dimension(:,:,:),pointer    :: p11     => null() ! (r,c,b) backscatter phase function, index 1
     real(kind=RKIND),dimension(:,:,:),pointer    :: p22     => null() ! (r,c,b) backscatter phase function, index 5

     integer,dimension(NRH_BINS):: rhi          ! pointer to rh LUT
     real(kind=RKIND),dimension(NRH_BINS):: rha ! slope on rh LUT

     contains
        procedure :: QueryByWavelength_1d
        procedure :: QueryByWavelength_2d
        procedure :: QueryByWavelength_3d
        procedure :: QueryByChannel_1d
        procedure :: QueryByChannel_2d
        procedure :: QueryByChannel_3d
        generic   :: Query => QueryByWavelength_1d, &
                              QueryByWavelength_2d, &
                              QueryByWavelength_3d, &
                              QueryByChannel_1d,    &
                              QueryByChannel_2d,    &
                              QueryByChannel_3d

        procedure:: getChannel
        procedure:: getWavelength
 end type GOCART2G_Mie

 interface GOCART2G_Mie
     module procedure GOCART2G_MieCreate
  end interface GOCART2G_Mie


 contains


!==================================================================================================================
 type(GOCART2G_Mie) function GOCART2G_MieCreate(dminfo,MieFile,wavelengths) result(self)
 use SMIOLf
#include "smiol_codes.inc"
!==================================================================================================================

!--- input arguments:
 type(dm_info),intent(in):: dminfo

 character(len=*),intent(in):: MieFile ! Mie table file name
 real(kind=RKIND),intent(in),dimension(:),optional:: wavelengths

!--- local variables needed to process the netCDF files:
 type(SMIOLf_context),pointer:: context
 type(SMIOLf_file),pointer   :: aop_file
 type(SMIOLf_decomp),pointer :: decomp   ! not used for non-decomposed variables

 logical:: l_gf,l_rhop
 integer:: i,imom,ip1,ipol,j,n,nn,nMom,nPol
 integer:: nb,nc,np,nr
 integer:: stat,ndims
 real(kind=RKIND):: yerr
 real(kind=RKIND),parameter:: undefval = 1.0e15
 real(kind=RKIND),dimension(:,:),pointer:: rhod,vol,area

!--- pointers used to read the netCDF variables using subroutines read_real_1d, read_real_2d, read_real_3d,
!    read_real_4d, and read_real_5d:
 integer(kind=I8KIND):: radius_size,rh_size,lambda_size,nMom_size,nPol_size

 real(kind=R4KIND),dimension(:),pointer        :: rh,lambda,radius,rLow,rUp
 real(kind=R4KIND),dimension(:,:),pointer      :: rEff,rMass,gf,rhop
 real(kind=R4KIND),dimension(:,:,:),pointer    :: qsca,qext,bsca,bext,g,bbck,refreal,refimag
 real(kind=R4KIND),dimension(:,:,:,:),pointer  :: pback

 !20240601_rkumar_start
 ! real(kind=R4KIND),dimension(:,:,:,:,:),pointer:: pmom
 !20240601_rkumar_end

!--- intermediate allocatable arrays used to interpolate the netCDF variables to the prescribed wavelengths:
 real(kind=RKIND),dimension(:),allocatable        :: lambda_r,input_r
 real(kind=RKIND),dimension(:,:,:),allocatable    :: bext_r,bsca_r,bbck_r,g_r,refreal_r,refimag_r
 real(kind=RKIND),dimension(:,:,:,:),allocatable  :: pback_r

 !20240601_rkumar_start
 ! real(kind=RKIND),dimension(:,:,:,:,:),allocatable:: pmom_r
 !20240601_rkumar_end

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter function GOCART2G_MieCreate:  '//MieFile)

 nullify(rh, lambda, radius, rLow, rUp)
 nullify(rEff, rMass, gf, rhop)
 nullify(qsca, qext, bsca, bext, g, bbck, refreal, refimag)
 nullify(pback)
 !20240601_rkumar_start
 ! nullify(pmom)
 !20240601_rkumar_end

!
!--- set up a context, within which SMIOL can do parallel I/O:
!
#ifdef MPAS_USE_MPI_F08
 stat = SMIOLf_init(dminfo%comm%mpi_val,1,dminfo%nprocs,context)
#else
 stat = SMIOLf_init(dminfo%comm,1,dminfo%nprocs,context)
#endif
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error initializing SMIOL context', messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)), messageType=MPAS_LOG_ERR)
    return
 endif


!
!--- open the netCDF file to be read:
!
 stat = SMIOLf_open_file(context,trim(MieFile),SMIOL_FILE_READ,aop_file)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error opening file '//trim(MieFile), messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)), messageType=MPAS_LOG_ERR)
    stat = SMIOLf_finalize(context)
    return
 endif


!
!--- inquire about the size of dimensions and initialize dimensions needed in GOCART2G_Mie:
!
 stat = SMIOLf_inquire_dim(aop_file,'radius',dimsize=radius_size)
 stat = SMIOLf_inquire_dim(aop_file,'rh',dimsize=rh_size)
 stat = SMIOLf_inquire_dim(aop_file,'lambda',dimsize=lambda_size)

 nMom = 0
 if(present(wavelengths)) then
    stat = SMIOLf_inquire_dim(aop_file,'nMom',dimsize=nMom_size)
    nMom = int(nMom_size)
 else
    nMom_size = 0
 endif

 nPol = 6
 if(present(wavelengths)) then
    stat = SMIOLf_inquire_dim(aop_file,'nPol',dimsize=nPol_size)
    nPol = int(nPol_size)
 else
    nPol_size = nPol
 endif

 self%nrh  = rh_size
 self%nbin = radius_size
 self%nMom = nMom
 self%nPol = nPol

 if(present(wavelengths)) then
    self%nch = size(wavelengths)
 else
    self%nch = lambda_size
 endif


!
!--- read variables as non-decomposed variables in the input netCDF files:
!    i.e., every MPI tasks reads the full variable:
!
 l_gf   = .false.
 l_rhop = .false.

 nullify(decomp)
 call read_real_1d(aop_file,decomp,'rh',rh)
 nullify(decomp)
 call read_real_1d(aop_file,decomp,'lambda',lambda)
 nullify(decomp)
 call read_real_1d(aop_file,decomp,'radius',radius)
 nullify(decomp)
 call read_real_1d(aop_file,decomp,'rLow',rLow)
 nullify(decomp)
 call read_real_1d(aop_file,decomp,'rUp',rUp)

 nullify(decomp)
 call read_real_2d(aop_file,decomp,'rEff',rEff)
 nullify(decomp)
 call read_real_2d(aop_file,decomp,'rMass',rMass)

 nullify(decomp)
 call read_real_3d(aop_file,decomp,'qsca',qsca)
 nullify(decomp)
 call read_real_3d(aop_file,decomp,'qext',qext)
 nullify(decomp)
 call read_real_3d(aop_file,decomp,'bsca',bsca)
 nullify(decomp)
 call read_real_3d(aop_file,decomp,'bext',bext)
 nullify(decomp)
 call read_real_3d(aop_file,decomp,'g',g)
 nullify(decomp)
 call read_real_3d(aop_file,decomp,'bbck',bbck)
 nullify(decomp)
 call read_real_3d(aop_file,decomp,'refreal',refreal)
 nullify(decomp)
 call read_real_3d(aop_file,decomp,'refimag',refimag)

 if(present(wavelengths)) then
    nullify(decomp)
    call read_real_4d(aop_file,decomp,'pback',pback)
    
    !20240601_rkumar_start
    ! nullify(decomp)
    ! call read_real_5d(aop_file,decomp,'pmom',pmom)
    !20240601_rkumar_end

    !growth factor:
    nullify(decomp)
    stat = SMIOLf_inquire_var(aop_file,'growth_factor',ndims=ndims)
    if(stat == SMIOL_SUCCESS) then
       l_gf = .true.
       call read_real_2d(aop_file,decomp,'growth_factor',gf)
    endif

    !wet particle density:
    nullify(decomp)
    stat = SMIOLf_inquire_var(aop_file,'rhop',ndims=ndims)
    if(stat == SMIOL_SUCCESS) then
       l_rhop = .true.
       call read_real_2d(aop_file,decomp,'rhop',rhop)
    endif
 endif

!--- close netCDF file:
 stat = SMIOLf_close_file(aop_file)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write(trim(SMIOLf_error_string(stat)), messageType=MPAS_LOG_ERR)
    stat = SMIOLf_finalize(context)
    return
 endif

 stat = SMIOLf_finalize(context)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error finalizing SMIOL context', messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)), messageType=MPAS_LOG_ERR)
    return
 endif


!
!--- initialize arrays contained in GOCART2G_Mie:
!
 if(.not.associated(self%wavelengths)) allocate(self%wavelengths(self%nch))
 if(present(wavelengths)) then
    self%wavelengths = wavelengths
 else
    self%wavelengths = lambda
 endif

 if(.not.associated(self%rh)      ) allocate(self%rh(self%nrh)                     )
 if(.not.associated(self%reff)    ) allocate(self%reff(self%nrh,self%nbin)         )
 if(.not.associated(self%bext)    ) allocate(self%bext(self%nrh,self%nch,self%nbin))
 if(.not.associated(self%bsca)    ) allocate(self%bsca(self%nrh,self%nch,self%nbin))
 if(.not.associated(self%bbck)    ) allocate(self%bbck(self%nrh,self%nch,self%nbin))
 if(.not.associated(self%g)       ) allocate(self%g(self%nrh,self%nch,self%nbin)   )
 if(.not.associated(self%refr)    ) allocate(self%refr(self%nrh,self%nch,self%nbin))
 if(.not.associated(self%refi)    ) allocate(self%refi(self%nrh,self%nch,self%nbin))
 if(.not.associated(self%p11)     ) allocate(self%p11(self%nrh,self%nch,self%nbin) )
 if(.not.associated(self%p22)     ) allocate(self%p22(self%nrh,self%nch,self%nbin) )

 if(present(wavelengths)) then
    if(.not.associated(self%gf)   ) allocate(self%gf(self%nrh,self%nbin)           )
    if(.not.associated(self%rhop) ) allocate(self%rhop(self%nrh,self%nbin)         )
    if(.not.associated(self%rhod) ) allocate(self%rhod(self%nrh,self%nbin)         )
    if(.not.associated(self%vol)  ) allocate(self%vol(self%nrh,self%nbin)          )
    if(.not.associated(self%area) ) allocate(self%area(self%nrh,self%nbin)         )

    if(.not.associated(self%pback)) allocate(self%pback(self%nrh,self%nch,self%nbin,self%nPol)         )
    
    !20240601_rkumar_start
    ! if(.not.associated(self%pmom) ) allocate(self%pmom(self%nrh,self%nch,self%nbin,self%nMom,self%nPol))
    !20240601_rkumar_end
 endif


!
!--- allocate intermediate arrays used to interpolate the netCDF variables to the prescribed wavelengths:
!
 if(.not.allocated(bext_r)    ) allocate(bext_r(int(lambda_size),self%nrh,self%nbin)   )
 if(.not.allocated(bsca_r)    ) allocate(bsca_r(int(lambda_size),self%nrh,self%nbin)   )
 if(.not.allocated(bbck_r)    ) allocate(bbck_r(int(lambda_size),self%nrh,self%nbin)   )
 if(.not.allocated(g_r)       ) allocate(g_r(int(lambda_size),self%nrh,self%nbin)      )
 if(.not.allocated(refreal_r) ) allocate(refreal_r(int(lambda_size),self%nrh,self%nbin))
 if(.not.allocated(refimag_r) ) allocate(refimag_r(int(lambda_size),self%nrh,self%nbin))
 if(present(wavelengths)) then
    if(.not.allocated(pback_r)) allocate(pback_r(int(lambda_size),self%nrh,self%nbin,self%nPol)         )
    
    !20240601_rkumar_start
    ! if(.not.allocated(pmom_r) ) allocate(pmom_r(int(lambda_size),self%nrh,self%nbin,self%nMom,self%nPol))
    !20240601_rkumar_end
 endif


!
!--- fill in the arrays contained in GOCART2G_Mie:
!
 self%rh   = real(rh,kind=RKIND)   ! relative humidity (fraction).
 self%reff = real(rEff,kind=RKIND) ! effective radius of bin (m).
 if(present(wavelengths)) then
    if(l_gf .or. l_rhop) then
       self%gf   = real(gf,kind=RKIND)
       self%rhop = real(rhop,kind=RKIND)
       self%rhod = undefval
       self%vol  = undefval
       self%area = undefval
    else
       self%gf   = undefval
       self%rhop = undefval
       self%rhod = undefval
       self%vol  = undefval
       self%area = undefval
    endif
 endif

 bext_r     = real(bext,kind=RKIND)
 bsca_r     = real(bsca,kind=RKIND)
 bbck_r     = real(bbck,kind=RKIND)
 g_r        = real(g,kind=RKIND)
 refreal_r  = real(refreal,kind=RKIND)
 refimag_r  = real(refimag,kind=RKIND)

 if (associated(bext))    then; deallocate(bext);    nullify(bext);    endif
 if (associated(bsca))    then; deallocate(bsca);    nullify(bsca);    endif
 if (associated(bbck))    then; deallocate(bbck);    nullify(bbck);    endif
 if (associated(g))       then; deallocate(g);       nullify(g);       endif
 if (associated(refreal)) then; deallocate(refreal); nullify(refreal); endif
 if (associated(refimag)) then; deallocate(refimag); nullify(refimag); endif

 if(present(wavelengths)) then
    pback_r = real(pback,kind=RKIND)
    
    !20240601_rkumar_start
    ! pmom_r  = real(pmom,kind=RKIND)
    !20240601_rkumar_end

    !20240601_rkumar_start
    if (associated(pback)) then; deallocate(pback); nullify(pback); endif
    ! if (associated(pmom))  then; deallocate(pmom);  nullify(pmom);  endif
    !20240601_rkumar_end
 endif


 if(present(wavelengths)) then
    if(.not.allocated(input_r) ) allocate(input_r(int(lambda_size)) )
    if(.not.allocated(lambda_r)) allocate(lambda_r(int(lambda_size)))
    do nn = 1,int(lambda_size)
       input_r(nn)  = 0._RKIND
       lambda_r(nn) = real(lambda(nn),kind=RKIND)
    enddo
    do j = 1,self%nbin
       do i = 1,self%nrh
          do n = 1,self%nch
             input_r(:) = bext_r(:,i,j)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%bext(i,n,j),yerr)
             input_r(:) = bsca_r(:,i,j)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%bsca(i,n,j),yerr)
             input_r(:) = bbck_r(:,i,j)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%bbck(i,n,j),yerr)
             input_r(:) = g_r(:,i,j)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%g(i,n,j)   ,yerr)
             input_r(:) = refreal_r(:,i,j)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%refr(i,n,j),yerr)
             input_r(:) = refimag_r(:,i,j)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%refi(i,n,j),yerr)

             do ipol = 1,self%nPol
                input_r(:) = pback_r(:,i,j,ipol)
                call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%pback(i,n,j,ipol),yerr)
             enddo

             !20240601_rkumar_start
             ! if(nMom > 0) then
             !    do imom = 1,self%nMom
             !       do ipol = 1,self%nPol
             !          input_r(:) = pmom_r(:,i,j,imom,ipol)
             !          call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n), &
             !                      self%pmom(i,n,j,imom,ipol),yerr)
             !       enddo
             !    enddo
             ! endif
             !20240601_rkumar_end

          enddo
       enddo
    enddo
    if(allocated(input_r) ) deallocate(input_r )
    if(allocated(lambda_r)) deallocate(lambda_r)
 else
    self%bext = reshape(bext_r   ,[self%nrh,int(lambda_size),self%nbin],order =[2,1,3])
    self%bsca = reshape(bsca_r   ,[self%nrh,int(lambda_size),self%nbin],order =[2,1,3])
    self%bbck = reshape(bbck_r   ,[self%nrh,int(lambda_size),self%nbin],order =[2,1,3])
    self%g    = reshape(g_r      ,[self%nrh,int(lambda_size),self%nbin],order =[2,1,3])
    self%refr = reshape(refreal_r,[self%nrh,int(lambda_size),self%nbin],order =[2,1,3])
    self%refi = reshape(refimag_r,[self%nrh,int(lambda_size),self%nbin],order =[2,1,3])
 endif

 if(allocated(bext_r)   ) deallocate(bext_r   )
 if(allocated(bsca_r)   ) deallocate(bsca_r   )
 if(allocated(bbck_r)   ) deallocate(bbck_r   )
 if(allocated(g_r)      ) deallocate(g_r      )
 if(allocated(refreal_r)) deallocate(refreal_r)
 if(allocated(refimag_r)) deallocate(refimag_r)
 if(present(wavelengths)) then
    if(allocated(pback_r)  ) deallocate(pback_r  )
    
    !20240601_rkumar_start
    ! if(allocated(pmom_r)   ) deallocate(pmom_r   )
    !20240601_rkumar_end
 endif

 if(present(wavelengths)) then
    self%p11 = self%pback(:,:,:,1)
    self%p22 = self%pback(:,:,:,5)
 else
    self%p11 = 1._RKIND
    self%p22 = 2._RKIND
 endif

 if(associated(rh))     then; deallocate(rh);     nullify(rh);     endif
 if(associated(lambda)) then; deallocate(lambda); nullify(lambda); endif
 if(associated(radius)) then; deallocate(radius); nullify(radius); endif
 if(associated(rLow))   then; deallocate(rLow);   nullify(rLow);   endif
 if(associated(rUp))    then; deallocate(rUp);    nullify(rUp);    endif
 if(associated(rEff))   then; deallocate(rEff);   nullify(rEff);   endif
 if(associated(rMass))  then; deallocate(rMass);  nullify(rMass);  endif
 if(associated(qsca))   then; deallocate(qsca);   nullify(qsca);   endif
 if(associated(qext))   then; deallocate(qext);   nullify(qext);   endif
 if(associated(gf))     then; deallocate(gf);     nullify(gf);     endif
 if(associated(rhop))   then; deallocate(rhop);   nullify(rhop);   endif

!
!--- remapping of RH to some high resolution representation for later use. RH input is scaled 0 - 0.99.
!    we resolve the map to 0 - 0.990 in steps of 0.001 (991 total steps):
 do j = 1, NRH_BINS
    do i = self%nrh, 1, -1
       if((j-1) .ge. int(self%rh(i)*1000)) then
          ip1 = i + 1
          self%rhi(j) = i
          if(ip1 .gt. self%nrh) then
             self%rha(j) = 0.
          else
             self%rha(j) = ( (j-1)/1000.-self%rh(i)) /  (self%rh(ip1)-self%rh(i))
          endif
          exit
       endif
    enddo
 enddo

 contains

   !-----------------------------------------------------------------------
   !  read_real_1d
   !-----------------------------------------------------------------------
   subroutine read_real_1d(file, decomp, varname, var)
      implicit none
      type (SMIOLf_file), intent(inout) :: file
      type (SMIOLf_decomp), pointer :: decomp
      character(len=*), intent(in) :: varname
      real(kind=R4KIND), dimension(:), pointer :: var

      double precision, dimension(:), pointer :: var_dbl
      integer :: vartype, ndims
      character(len=256), dimension(1) :: dimname
      integer(kind=SMIOL_offset_kind) :: dimsize
      integer :: ierr

      nullify(var)
      ierr = SMIOLf_inquire_var(file, varname, vartype=vartype, ndims=ndims)
      if (ierr /= SMIOL_SUCCESS .or. ndims /= 1) return

      ierr = SMIOLf_inquire_var(file, varname, dimnames=dimname)
      ierr = SMIOLf_inquire_dim(file, dimname(1), dimsize=dimsize)

      if (vartype == SMIOL_REAL32) then
          allocate(var(dimsize))
          ierr = SMIOLf_get_var(file, varname, decomp, var)
      else if (vartype == SMIOL_REAL64) then
          allocate(var(dimsize))
          allocate(var_dbl(dimsize))
          ierr = SMIOLf_get_var(file, varname, decomp, var_dbl)
          var(:) = var_dbl(:)
          deallocate(var_dbl)
      end if
   end subroutine read_real_1d

   !-----------------------------------------------------------------------
   !  read_real_2d
   !-----------------------------------------------------------------------
   subroutine read_real_2d(file, decomp, varname, var)
      implicit none
      type (SMIOLf_file), intent(inout) :: file
      type (SMIOLf_decomp), pointer :: decomp
      character(len=*), intent(in) :: varname
      real(kind=R4KIND), dimension(:,:), pointer :: var

      double precision, dimension(:,:), pointer :: var_dbl
      integer :: n, vartype, ndims
      character(len=256), dimension(2) :: dimname
      integer(kind=SMIOL_offset_kind), dimension(2) :: dimsize
      integer :: ierr

      nullify(var)
      ierr = SMIOLf_inquire_var(file, varname, vartype=vartype, ndims=ndims)
      if (ierr /= SMIOL_SUCCESS .or. ndims /= 2) return

      ierr = SMIOLf_inquire_var(file, varname, dimnames=dimname)
      do n = 1, 2
          ierr = SMIOLf_inquire_dim(file, dimname(n), dimsize=dimsize(n))
      enddo

      if (vartype == SMIOL_REAL32) then
          allocate(var(dimsize(1),dimsize(2)))
          ierr = SMIOLf_get_var(file, varname, decomp, var)
      else if (vartype == SMIOL_REAL64) then
          allocate(var(dimsize(1),dimsize(2)))
          allocate(var_dbl(dimsize(1),dimsize(2)))
          ierr = SMIOLf_get_var(file, varname, decomp, var_dbl)
          var(:,:) = var_dbl(:,:)
          deallocate(var_dbl)
      end if
   end subroutine read_real_2d

   !-----------------------------------------------------------------------
   !  read_real_3d
   !-----------------------------------------------------------------------
   subroutine read_real_3d(file, decomp, varname, var)
      implicit none
      type (SMIOLf_file), intent(inout) :: file
      type (SMIOLf_decomp), pointer :: decomp
      character(len=*), intent(in) :: varname
      real(kind=R4KIND), dimension(:,:,:), pointer :: var

      double precision, dimension(:,:,:), pointer :: var_dbl
      integer :: n, vartype, ndims
      character(len=256), dimension(3) :: dimname
      integer(kind=SMIOL_offset_kind), dimension(3) :: dimsize
      integer :: ierr

      nullify(var)
      ierr = SMIOLf_inquire_var(file, varname, vartype=vartype, ndims=ndims)
      if (ierr /= SMIOL_SUCCESS .or. ndims /= 3) return

      ierr = SMIOLf_inquire_var(file, varname, dimnames=dimname)
      do n = 1, 3
          ierr = SMIOLf_inquire_dim(file, dimname(n), dimsize=dimsize(n))
      enddo

      if (vartype == SMIOL_REAL32) then
          allocate(var(dimsize(1),dimsize(2),dimsize(3)))
          ierr = SMIOLf_get_var(file, varname, decomp, var)
      else if (vartype == SMIOL_REAL64) then
          allocate(var(dimsize(1),dimsize(2),dimsize(3)))
          allocate(var_dbl(dimsize(1),dimsize(2),dimsize(3)))
          ierr = SMIOLf_get_var(file, varname, decomp, var_dbl)
          var(:,:,:) = var_dbl(:,:,:)
          deallocate(var_dbl)
      end if
   end subroutine read_real_3d

   !-----------------------------------------------------------------------
   !  read_real_4d
   !-----------------------------------------------------------------------
   subroutine read_real_4d(file, decomp, varname, var)
      implicit none
      type (SMIOLf_file), intent(inout) :: file
      type (SMIOLf_decomp), pointer :: decomp
      character(len=*), intent(in) :: varname
      real(kind=R4KIND), dimension(:,:,:,:), pointer :: var

      double precision, dimension(:,:,:,:), pointer :: var_dbl
      integer :: n, vartype, ndims
      character(len=256), dimension(4) :: dimname
      integer(kind=SMIOL_offset_kind), dimension(4) :: dimsize
      integer :: ierr

      nullify(var)
      ierr = SMIOLf_inquire_var(file, varname, vartype=vartype, ndims=ndims)
      if (ierr /= SMIOL_SUCCESS .or. ndims /= 4) return

      ierr = SMIOLf_inquire_var(file, varname, dimnames=dimname)
      do n = 1, 4
          ierr = SMIOLf_inquire_dim(file, dimname(n), dimsize=dimsize(n))
      enddo

      if (vartype == SMIOL_REAL32) then
          allocate(var(dimsize(1),dimsize(2),dimsize(3),dimsize(4)))
          ierr = SMIOLf_get_var(file, varname, decomp, var)
      else if (vartype == SMIOL_REAL64) then
          allocate(var(dimsize(1),dimsize(2),dimsize(3),dimsize(4)))
          allocate(var_dbl(dimsize(1),dimsize(2),dimsize(3),dimsize(4)))
          ierr = SMIOLf_get_var(file, varname, decomp, var_dbl)
          var(:,:,:,:) = var_dbl(:,:,:,:)
          deallocate(var_dbl)
      end if
   end subroutine read_real_4d

   !-----------------------------------------------------------------------
   !  read_real_5d
   !-----------------------------------------------------------------------
   subroutine read_real_5d(file, decomp, varname, var)
      implicit none
      type (SMIOLf_file), intent(inout) :: file
      type (SMIOLf_decomp), pointer :: decomp
      character(len=*), intent(in) :: varname
      real(kind=R4KIND), dimension(:,:,:,:,:), pointer :: var

      double precision, dimension(:,:,:,:,:), pointer :: var_dbl
      integer :: n, vartype, ndims
      character(len=256), dimension(5) :: dimname
      integer(kind=SMIOL_offset_kind), dimension(5) :: dimsize
      integer :: ierr

      nullify(var)
      ierr = SMIOLf_inquire_var(file, varname, vartype=vartype, ndims=ndims)
      if (ierr /= SMIOL_SUCCESS .or. ndims /= 5) return

      ierr = SMIOLf_inquire_var(file, varname, dimnames=dimname)
      do n = 1, 5
          ierr = SMIOLf_inquire_dim(file, dimname(n), dimsize=dimsize(n))
      enddo

      if (vartype == SMIOL_REAL32) then
          allocate(var(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5)))
          ierr = SMIOLf_get_var(file, varname, decomp, var)
      else if (vartype == SMIOL_REAL64) then
          allocate(var(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5)))
          allocate(var_dbl(dimsize(1),dimsize(2),dimsize(3),dimsize(4),dimsize(5)))
          ierr = SMIOLf_get_var(file, varname, decomp, var_dbl)
          var(:,:,:,:,:) = var_dbl(:,:,:,:,:)
          deallocate(var_dbl)
      end if
   end subroutine read_real_5d

 end function GOCART2G_MieCreate

!==================================================================================================================
 subroutine polint(x,y,n,xWant,yWant,yErr)
 implicit none
 integer,intent(in):: n
 real(kind=RKIND),intent(in):: x(n),y(n)
 real(kind=RKIND),intent(inout):: xWant,yWant,yErr

 character(len=255):: msg
 integer:: i, j
 real(kind=RKIND):: dx, slope

 i = 0
 if(xWant .lt. x(1)) i = 1
 if(xWant .gt. x(n)) i = n

 if(i .eq. 0) then
    do j = 1, n
       if(xWant .ge. x(j)) i = j
    enddo
 endif

 if(i .eq. n) then
    slope = 0.
 else
    slope = (y(i+1)-y(i)) / (x(i+1)-x(i))
 endif
 dx  = xWant - x(i)
 yWant = y(i) + slope*dx
 yErr = 0.

 end subroutine polint

!==================================================================================================================
#define RANK_ 1
#include "MieQuery.H"
#undef RANK_

#define RANK_ 2
#include "MieQuery.H"
#undef RANK_

#define RANK_ 3
#include "MieQuery.H"
#undef RANK_

!==================================================================================================================
 integer function getChannel(this,wavelength,rc) result (ch)
 class(GOCART2G_Mie),intent(in):: this
 real(kind=RKIND),intent(in):: wavelength
 integer,intent(out),optional:: rc

 integer:: i
 real(kind=RKIND),parameter:: w_tol = 1.e-9

 ch = -1
 do i = 1, this%nch
    if (abs(this%wavelengths(i)-wavelength) <= w_tol) then
        ch = i
        exit
    endif
 enddo

 if (present(rc)) rc = 0

 if (ch < 0) then
    !$omp critical (GetCha)
    print*, "wavelength of ",wavelength, " is an invalid value."
    !$omp end critical (GetCha)
    if (present(rc)) rc = -1
 endif

 end function getChannel

!==================================================================================================================
 real function getWavelength(this,ith_channel,rc) result (wavelength)
 class(GOCART2G_Mie),intent(in):: this
 integer,intent(in):: ith_channel
 integer,intent(out),optional:: rc

 if (present(rc)) rc = 0

 if (ith_channel <=0 .or. ith_channel > this%nch ) then
    !$omp critical (GetWav)
    print*, "The channel of ",ith_channel, " is an invalid channel number."
    !$omp end critical (GetWav)
    if (present(rc)) rc = -1
    wavelength = -1. 
    return
 endif

 wavelength = this%wavelengths(ith_channel)

 end function getWavelength

 end module gocart2G_MieMod_smiol
