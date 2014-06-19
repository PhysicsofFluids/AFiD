      subroutine openfi
      use param
      implicit none

!EP    nusse.out  in vmaxv.F
      open(95,file='nusse.out',status='unknown',access='sequential',
     % position='append')

!EP    nusse2.out  in densmc.F
       open(97,file="nusse2.out",status='unknown',access='sequential',
     % position='append')

!EP   rms_vel.out in stst.F
       open(94,file='rms_vel.out',status='unknown',position='append',
     & access='sequential')

#ifdef BALANCE
!EP   nusse3.out in balance.F
      open(92,file='nusse3.out',status='unknown',access='sequential',
     & position='append')
#endif

#ifdef STATS2
!EP   probeline1.out in stst2.F
       open(44,file='probeline1.out',status='unknown',position='append',
     & access='sequential')
      
!EP   probeline2.out in stst2.F
       open(45,file='probeline2.out',status='unknown',position='append',
     & access='sequential')

!EP   ring1.out in stst2.F
       open(46,file='ring1.out',status='unknown',position='append',
     & access='sequential')

!EP   ring2.out in stst2.F
       open(47,file='ring2.out',status='unknown',position='append',
     & access='sequential')
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
#ifdef STATS2
      close(44)
      close(45)
      close(46)
      close(47)
#endif
      return 
      end
