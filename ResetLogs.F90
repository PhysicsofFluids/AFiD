!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ResetLogs.F90                                  !
!    CONTAINS: subroutine ResetLogs                       !
!                                                         ! 
!    PURPOSE: Initialization routine. Reset all log files !
!     if necessary                                        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ResetLogs
      use param
      implicit none

      if(resetlogstime) then    

!EP    nusse.out  in GlobalQuantities.F
      open(95,file='nu_vol.out',status='unknown',access='sequential', &
        position='append')
      close(95,status='delete')

!EP    nusse2.out  in CalculatePlateNu.F
       open(97,file="nu_plate.out",status='unknown',access='sequential', &
        position='append')
      close(97,status='delete')

!EP   nusse3.out in CalcDissipationNu.F
      if (disscal) then
      open(92,file='nu_diss.out',status='unknown',access='sequential', &
       position='append')
      close(92,status='delete')
      end if

!EP   rms_vel.out in GlobalQuantities.F
       open(94,file='rms_vel.out',status='unknown',position='append', &
        access='sequential')
      close(94,status='delete')


      endif

      return
      end   
      

