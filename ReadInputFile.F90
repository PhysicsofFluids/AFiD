!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ReadInputFile.F90                              !
!    CONTAINS: subroutine ReadInputFile                   !
!                                                         ! 
!    PURPOSE: Read parameters from bou.in file            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ReadInputFile
      use param
      implicit none
      character(len=4) :: dummy
      integer flagstat,flagbal,stst3flag
      logical fexist

      open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) nxm,nym,nzm,nsst,nread
        read(15,301) dummy
        read(15,*) ntst,walltimemax,tout,tmax,ireset
        read(15,301) dummy
        read(15,*) alx3,istr3,str3
        read(15,301) dummy
        read(15,*) ylen,zlen
        read(15,301) dummy
        read(15,*) ray,pra,dt,resid,limitCFL
        read(15,301) dummy
        read(15,*) flagstat,flagbal,tsta,starea
        read(15,301) dummy
        read(15,*) inslws,inslwn
        read(15,301) dummy       
        read(15,*) idtv,dtmin,dtmax,limitVel
        read(15,301) dummy       
        read(15,*) tframe,stst3flag
301     format(a4)                
      close(15)

      nx=nxm+1
      ny=nym+1                                                          
      nz=nzm+1                                                          

!m============================================
!
!     DEFINITIONS FOR THE NATURAL CONVECTION
!
      ren = dsqrt(ray/pra)
      pec = dsqrt(pra*ray)
      pi=2.d0*dasin(1.d0)                          
!                                                                       
!
      if(flagstat.ne.0) statcal = .true.
      if(idtv.eq.0) variabletstep = .false.
      if(flagbal.ne.0) disscal = .true.
      if(nread.ne.0) readflow = .true.
      if(ireset.ne.0) resetlogstime = .true.

      if(stst3flag.ne.0) then
       inquire(file='./stst3.in', exist=fexist) 
       if(fexist) then
        dumpslabs = .true.
       else
        write(6,*) "stst3.in not found, turning off slab dump"
       end if
      endif


      if(starea.ne.0) then 
        readstats = .true.
        if (.not. readflow) write(6,*) 'Warning: Restarting flowfield with statistics read'
      endif


      return 
      end
