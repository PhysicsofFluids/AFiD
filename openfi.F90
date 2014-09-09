!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: openfi.F90                                     !
!    CONTAINS: subroutine openfi,closefi                  !
!                                                         ! 
!    PURPOSE: Initialization routine. Open/close all log  !
!     files and rewind if necessary                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine openfi
      use param
      implicit none

!EP    nusse.out  in globalquantities.F
      open(95,file='nusse.out',status='unknown',access='sequential', &
        position='append')

!EP    nusse2.out  in densmc.F
       open(97,file="nusse2.out",status='unknown',access='sequential', &
        position='append')

!EP   rms_vel.out in globalquantities.F
       open(94,file='rms_vel.out',status='unknown',position='append', &
        access='sequential')

!EP   nusse3.out in balance.F
      if (balcal) then
      open(92,file='nusse3.out',status='unknown',access='sequential', &
       position='append')
      end if

      if(ireset.eq.1) then    
      rewind(92)
      rewind(94)
      rewind(95)
      rewind(97)
      endif

      return
      end   
      


      subroutine closefi
      use param
      implicit none
      close(95)
      close(97)
      close(94)
      if (balcal) close(92)
      return 
      end
