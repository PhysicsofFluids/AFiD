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

#ifdef BALANCE
!EP   nusse3.out in balance.F
      open(92,file='nusse3.out',status='unknown',access='sequential', &
       position='append')
#endif

      if(ireset.eq.1) then    
      rewind(92)
      rewind(94)
      rewind(95)
      rewind(97)
      endif

      return
      end   
      


      subroutine closefi
      implicit none
      close(95)
      close(97)
      close(94)
#ifdef BALANCE
      close(92)
#endif
      return 
      end
