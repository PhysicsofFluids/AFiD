      subroutine ReadInputFile
      use param
      implicit none
      character(len=4) :: dummy
      integer flagstat,flagbal,stst3flag

      open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) n1,n2,n3,nsst,nwrit,nread
        read(15,301) dummy
        read(15,*) ntst,walltimemax,tpin,tmax,ireset
        read(15,301) dummy
        read(15,*) alx3,istr3,str3
        read(15,301) dummy
        read(15,*) rext,rext2
        read(15,301) dummy
        read(15,*) ray,pra,dt,resid,cflmax
        read(15,301) dummy
        read(15,*) flagstat,flagbal,tsta,starea
        read(15,301) dummy
        read(15,*) inslws,inslwn
        read(15,301) dummy       
        read(15,*) idtv,dtmin,dtmax,cfllim,vlim
        read(15,301) dummy       
        read(15,*) tframe,stst3flag
301     format(a4)                
      close(15)

      n1m=n1-1                                                          
      n2m=n2-1                                                          
      n3m=n3-1

!m============================================
!
!     DEFINITIONS FOR THE NATURAL CONVECTION
!
      ren = dsqrt(ray/pra)
      pec = dsqrt(pra*ray)
      pi=2.d0*dasin(1.d0)                          
!                                                                       
!
      if(flagstat.eq.0) then
        statcal = .false.
      else
        statcal = .true.
      endif

      if(flagbal.eq.0) then
        balcal = .false.
      else
        balcal = .true.
      endif

      if(stst3flag.eq.0) then
        dumpslabs = .false.
      else
        dumpslabs = .true.
      endif

      return 
      end
