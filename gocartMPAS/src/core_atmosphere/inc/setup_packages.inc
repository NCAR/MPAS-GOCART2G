#ifdef MPAS_DEBUG
#define COMMA ,
#define PACKAGE_LOGIC_PRINT(M) call mpas_log_write(M)
#else
#define PACKAGE_LOGIC_PRINT(M) ! M
#endif

   !
   ! WARNING: This function is automatically generated at compile time.
   !          Any modifications to this code will be lost when MPAS is recompiled.
   !
   function atm_setup_packages_when(configPool, packagePool) result(ierr)

      use mpas_derived_types, only : mpas_pool_type
      use mpas_log, only : mpas_log_write

      implicit none

      integer :: ierr

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool


      ierr = 0

      call mpas_log_write('')
      call mpas_log_write('Configuring registry-specified packages...')

      call atm_setup_les_package(configPool, packagePool)
      call atm_setup_sfclayer_package(configPool, packagePool)

      call mpas_log_write('----- done configuring registry-specified packages -----')
      call mpas_log_write('')

   end function atm_setup_packages_when


   !
   ! WARNING: This subroutine is automatically generated at compile time.
   !          Any modifications to this code will be lost when MPAS is recompiled.
   !
   subroutine atm_setup_les_package(configPool, packagePool)

      use mpas_kind_types, only : RKIND, StrKIND
      use mpas_derived_types, only : mpas_pool_type
      use mpas_log, only : mpas_log_write

      implicit none

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool

      logical, pointer :: lesActive

      character(len=StrKIND), pointer :: config_les_model

      nullify(config_les_model)
      call mpas_pool_get_config(configPool, 'config_les_model', config_les_model)

      nullify(lesActive)
      call mpas_pool_get_package(packagePool, 'lesActive', lesActive)

      PACKAGE_LOGIC_PRINT('')
      PACKAGE_LOGIC_PRINT("  les is active when (trim(config_les_model) /= 'none')")
      PACKAGE_LOGIC_PRINT('    namelist settings:')
      PACKAGE_LOGIC_PRINT('    ------------------')
      PACKAGE_LOGIC_PRINT('    config_les_model = '//trim(config_les_model))

      lesActive = ( trim(config_les_model) /= 'none' )
      call mpas_log_write('  les : $l', logicArgs=[lesActive])

   end subroutine atm_setup_les_package


   !
   ! WARNING: This subroutine is automatically generated at compile time.
   !          Any modifications to this code will be lost when MPAS is recompiled.
   !
   subroutine atm_setup_sfclayer_package(configPool, packagePool)

      use mpas_kind_types, only : RKIND, StrKIND
      use mpas_derived_types, only : mpas_pool_type
      use mpas_log, only : mpas_log_write

      implicit none

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool

      logical, pointer :: sfclayerActive

      character(len=StrKIND), pointer :: config_physics_suite
      character(len=StrKIND), pointer :: config_sfclayer_scheme

      nullify(config_physics_suite)
      call mpas_pool_get_config(configPool, 'config_physics_suite', config_physics_suite)
      nullify(config_sfclayer_scheme)
      call mpas_pool_get_config(configPool, 'config_sfclayer_scheme', config_sfclayer_scheme)

      nullify(sfclayerActive)
      call mpas_pool_get_package(packagePool, 'sfclayerActive', sfclayerActive)

      PACKAGE_LOGIC_PRINT('')
      PACKAGE_LOGIC_PRINT("  sfclayer is active when ((trim(config_physics_suite) /= 'none' .and. trim(config_sfclayer_scheme) == &
&'suite') .or.  any(trim(config_sfclayer_scheme) == [ character(len=StrKind) :: 'sf_monin_obukhov', 'sf_monin_obukhov_rev', &
&'sf_mynn' ]))")
      PACKAGE_LOGIC_PRINT('    namelist settings:')
      PACKAGE_LOGIC_PRINT('    ------------------')
      PACKAGE_LOGIC_PRINT('    config_physics_suite = '//trim(config_physics_suite))
      PACKAGE_LOGIC_PRINT('    config_sfclayer_scheme = '//trim(config_sfclayer_scheme))

      sfclayerActive = ( (trim(config_physics_suite) /= 'none' .and. trim(config_sfclayer_scheme) == 'suite') .or.  &
any(trim(config_sfclayer_scheme) == [ character(len=StrKind) :: 'sf_monin_obukhov', 'sf_monin_obukhov_rev', 'sf_mynn' ]) )
      call mpas_log_write('  sfclayer : $l', logicArgs=[sfclayerActive])

   end subroutine atm_setup_sfclayer_package

