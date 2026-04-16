! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
!Initialization of alarms for PRM and PRM input data. 
!Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-04-02.
!
!==================================================================================================================
 module mpas_chemistry_plumerisemodel_manager

 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines,only: mpas_pool_get_config,mpas_pool_get_subpool
 use mpas_stream_manager
 use mpas_timekeeping


 implicit none
    
    private
    public:: plumerisemodel_timetracker_init

    !--- Global variables

        integer,public:: iyear
        integer,public:: imonth
        integer,public:: iday
        integer,public:: ihour
        integer,public:: iminute
        integer,public:: isecond

        character(len=*),parameter:: plumerisemodel_lbc_AlarmID = 'plumerisemodel_lbc'
        character(len=*),parameter:: plumerisemodel_AlarmID     = 'plumerisemodel'

 contains

!==================================================================================================================
 subroutine plumerisemodel_timetracker_init(configs, clock)
!==================================================================================================================


    !*******************************************************************************
    !implicit none
    !*******************************************************************************

    !--- input arguments:
        type(mpas_pool_type),intent(in):: configs


    !--- local variables:

        type(MPAS_Clock_type),intent(inout):: clock
        type(MPAS_Time_Type):: startTime, alarmStartTime
        type(MPAS_TimeInterval_Type):: alarmTimeStep

        character(len=StrKIND),pointer:: prm_lbc_interval
        character(len=StrKIND),pointer:: prm_interval

        !Debug Tracking Variables
        !character(len=StrKIND),pointer:: printalarmStartTime

        integer:: ierr


    call mpas_log_write(' ')
    call mpas_log_write('--- enter subroutine plumerisemodel_timetracker_init:')


    !*******************************************************************************
    ! Computes startTime
    !*******************************************************************************

        !--- Get start time
                
            ierr = 0
            startTime = mpas_get_clock_time(clock,MPAS_START_TIME,ierr)
            if(ierr /=0) &
                call mpas_log_write('--- plumerisemodel_timetracker_init: error getting startTime', &
                                        messageType=MPAS_LOG_CRIT)

        !--- Retrieve year, month, etc
            ierr = 0
            call mpas_get_time(curr_time=startTime,YYYY=iyear,MM=imonth,DD=iday,H=ihour,M=iminute,S=isecond,ierr=ierr)
            if(ierr /=0) then
                call mpas_log_write('--- plumerisemodel_timetracker_init: error deriving iyear,imonth,etc...', &
                                    messageType=MPAS_LOG_CRIT)
            else
                call mpas_log_write('    startTime:')
                call mpas_log_write('    YEAR   = $i',intArgs=(/iyear/))
                call mpas_log_write('    MONTH  = $i',intArgs=(/imonth/))
                call mpas_log_write('    DAY    = $i',intArgs=(/iday/))
                call mpas_log_write('    HOUR   = $i',intArgs=(/ihour/))
                call mpas_log_write('    MINUTE = $i',intArgs=(/iminute/))
                call mpas_log_write('    SECOND = $i',intArgs=(/isecond/))
            endif


    !*******************************************************************************
    ! Initializes alarms to update PRM :
    !*******************************************************************************

        !---  Initialize alarms to update PRM Input data (Fire Size and FRP):

            call mpas_pool_get_config(configs,'config_prm_lbc_interval',prm_lbc_interval)
            ierr = 0
            call mpas_set_timeInterval(alarmTimeStep,timeString=prm_lbc_interval,ierr=ierr)
            alarmStartTime = startTime
            call mpas_add_clock_alarm(clock,plumerisemodel_lbc_AlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
            if(ierr /= 0) &
            call mpas_log_write('--- plumerisemodel_timetracker_init: error creating plumerisemodel_lbc_AlarmID', &
                                    messageType=MPAS_LOG_CRIT)


            !--- Debug Tracking Variables
            call mpas_log_write('    prm_lbc_interval is  '//trim(prm_lbc_interval))
            
        
        !---  Initialize alarms to update PRM :

            call mpas_pool_get_config(configs,'config_prm_interval',prm_interval)
            ierr = 0
            call mpas_set_timeInterval(alarmTimeStep,timeString=prm_interval,ierr=ierr)
            alarmStartTime = startTime
            call mpas_add_clock_alarm(clock,plumerisemodel_AlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
            if(ierr /= 0) &
            call mpas_log_write('--- plumerisemodel_timetracker_init: error creating plumerisemodel_AlarmID', &
                                    messageType=MPAS_LOG_CRIT)


            !--- Debug Tracking Variables
            call mpas_log_write('    prm_interval is  '//trim(prm_interval))


 call mpas_log_write('--- end subroutine plumerisemodel_timetracker_init.')
 call mpas_log_write(' ')
 
 end subroutine plumerisemodel_timetracker_init
!==================================================================================================================




!==================================================================================================================
 end module mpas_chemistry_plumerisemodel_manager
!==================================================================================================================
