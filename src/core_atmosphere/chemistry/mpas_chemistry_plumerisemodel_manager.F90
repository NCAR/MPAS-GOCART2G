! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!================================================================================================================== 
!
!Initialization of alarms for PRM and PRM input data. 
!Update of PRM input data (Fire Size, FRP) 
!Run PRM to calculate smoke injection height
!Jaqueline Pereira (jaqueline.pereira@inpe.br) / 2026-04-02.
!
!==================================================================================================================
 module mpas_chemistry_plumerisemodel_manager

 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines,only: mpas_pool_get_config,mpas_pool_get_subpool
 use mpas_stream_manager
 use mpas_timekeeping

 use mpas_chemistry_plumerisemodel_lbc_update
 use mpas_chemistry_plumerisemodel_lbc_init
 use mpas_chemistry_plumerisemodel_driver


 implicit none
    
    private
    public:: plumerisemodel_timetracker_init, plumerisemodel_timetracker_lbc, plumerisemodel_timetracker_prm

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
 subroutine plumerisemodel_timetracker_lbc(domain,clock,stream_manager)
!==================================================================================================================


    !*******************************************************************************
    !implicit none
    !*******************************************************************************

        !--- inout arguments:
            
            type(MPAS_Clock_type),intent(inout):: clock
            type(domain_type),intent(inout):: domain
            type(MPAS_streamManager_type),intent(inout):: stream_manager
                    

        !--- local variables:
            
            type(block_type),pointer:: block
            type(MPAS_Time_Type):: currTime
            type(mpas_pool_type),pointer:: mesh
            type(mpas_pool_type),pointer:: FINN_biob_prm_lbc
            type(mpas_pool_type),pointer:: biob_prm_lbc
            
            character(len=StrKIND):: timeStamp
            character(len=StrKIND),pointer:: prm_lbc_interval

            integer:: ierr

            

        call mpas_log_write(' ')
        call mpas_log_write('--- enter subroutine plumerisemodel_timetracker_lbc:')
        call mpas_log_write(' ')

    !*******************************************************************************
    ! Computes Current Time
    !*******************************************************************************

        !--- Get current time
                
            currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)


        !--- Retrieve current year, month, etc

            call mpas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=iyear,MM=imonth, &
                        DD=iday,H=ihour,M=iminute,S=isecond,ierr=ierr)
            
            call mpas_log_write('    currTime:')
            call mpas_log_write('    YEAR   = $i',intArgs=(/iyear/))
            call mpas_log_write('    MONTH  = $i',intArgs=(/imonth/))
            call mpas_log_write('    DAY    = $i',intArgs=(/iday/))
            call mpas_log_write('    HOUR   = $i',intArgs=(/ihour/))
            call mpas_log_write('    MINUTE = $i',intArgs=(/iminute/))
            call mpas_log_write('    SECOND = $i',intArgs=(/isecond/))


    !*******************************************************************************
    ! Update PRM input data (Fire Size, FRP):
    !*******************************************************************************

        !--- check to see if it is time to update PRM input data:
            
            block => domain%blocklist
            do while(associated(block))

                call mpas_pool_get_config(domain%blocklist%configs,'config_prm_lbc_interval',prm_lbc_interval)
                call mpas_log_write('    config_prm_lbc_interval = '//prm_lbc_interval)
    

                call mpas_pool_get_subpool(block%structs,'mesh'                ,mesh              )
                call mpas_pool_get_subpool(block%structs,'FINN_biob_prm_lbc'   ,FINN_biob_prm_lbc )
                call mpas_pool_get_subpool(block%structs,'biob_prm_lbc'        ,biob_prm_lbc      )
                

                if(mpas_is_alarm_ringing(clock,plumerisemodel_lbc_AlarmID,ierr=ierr)) then
                    call mpas_reset_clock_alarm(clock, plumerisemodel_lbc_AlarmID, ierr=ierr)
                    call mpas_log_write('    time to update PRM Input Data:')
                    call plumerisemodel_lbc_update(clock,stream_manager,mesh,FINN_biob_prm_lbc, biob_prm_lbc)
                else
                    call mpas_log_write('    No need to update PRM Input Data')
                endif

                !--- Fills local arrays with updated values:
                call plumerisemodel_lbc_init(mesh, biob_prm_lbc)

            block => block%next
            end do


call mpas_log_write('--- end subroutine plumerisemodel_timetracker_lbc.')
call mpas_log_write(' ')

end subroutine plumerisemodel_timetracker_lbc
!==================================================================================================================


!==================================================================================================================
 subroutine plumerisemodel_timetracker_prm(domain,clock)
!==================================================================================================================


    !*******************************************************************************
    !implicit none
    !*******************************************************************************

        !--- inout arguments:
            
            type(MPAS_Clock_type),intent(inout):: clock
            type(domain_type),intent(inout):: domain
                    

        !--- local variables:

            type(block_type),pointer:: block
            type(MPAS_Time_Type):: currTime
            type(mpas_pool_type),pointer:: mesh

            character(len=StrKIND):: timeStamp
            character(len=StrKIND),pointer:: prm_interval

            integer:: ierr

        call mpas_log_write(' ')
        call mpas_log_write('--- enter subroutine plumerisemodel_timetracker_prm:')



    !*******************************************************************************
    ! Computes Current Time
    !*******************************************************************************

        !--- Get current time
                
            currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)


        !--- Retrieve current year, month, etc

            call mpas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=iyear,MM=imonth, &
                        DD=iday,H=ihour,M=iminute,S=isecond,ierr=ierr)
            
            call mpas_log_write('    currTime:')
            call mpas_log_write('    YEAR   = $i',intArgs=(/iyear/))
            call mpas_log_write('    MONTH  = $i',intArgs=(/imonth/))
            call mpas_log_write('    DAY    = $i',intArgs=(/iday/))
            call mpas_log_write('    HOUR   = $i',intArgs=(/ihour/))
            call mpas_log_write('    MINUTE = $i',intArgs=(/iminute/))
            call mpas_log_write('    SECOND = $i',intArgs=(/isecond/))

    !*******************************************************************************
    ! Run PRM:
    !*******************************************************************************

        !--- check to see if it is time to run PRM:
            
            block => domain%blocklist
            do while(associated(block))

                call mpas_pool_get_config(domain%blocklist%configs,'config_prm_interval',prm_interval)
                call mpas_log_write('    config_prm_interval = '//prm_interval)
    

                call mpas_pool_get_subpool(block%structs,'mesh'                ,mesh              )              

                if(mpas_is_alarm_ringing(clock,plumerisemodel_AlarmID,ierr=ierr)) then
                    call mpas_reset_clock_alarm(clock, plumerisemodel_AlarmID, ierr=ierr)
                    call mpas_log_write('    time to run PRM:')
                    call plumerisemodel_driver(domain, mesh)
                else
                    call mpas_log_write('    No need to run PRM')
                endif

            block => block%next
            end do

call mpas_log_write('--- end subroutine plumerisemodel_timetracker_prm.')
call mpas_log_write(' ')

end subroutine plumerisemodel_timetracker_prm
!==================================================================================================================



!==================================================================================================================
 end module mpas_chemistry_plumerisemodel_manager
!==================================================================================================================
