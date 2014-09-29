!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitializeTimeMarchScheme.F90                  !
!    CONTAINS: subroutine InitializeTimeMarchScheme       !
!                                                         ! 
!    PURPOSE: Initialize the time-marching constants for  !
!     the integrator                                      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine InitializeTimeMarchScheme
      use param
      implicit none
      integer ns

      if(nsst.gt.1) then   
        gam(1)=8.d0/15.d0
        gam(2)=5.d0/12.d0
        gam(3)=3.d0/4.d0
        rom(1)=0.d0
        rom(2)=-17.d0/60.d0
        rom(3)=-5.d0/12.d0
!m======================================================
      if(ismaster) then
        write(6,100) (gam(ns),ns=1,nsst),(rom(ns),ns=1,nsst)
  100   format(/,5x,'The time scheme is a III order Runge-Kutta' &
        ,4x,'gam= ',3f8.3,4x,'ro= ',3f8.3)
      endif
!m======================================================
      else                                                              
        gam(1)=1.5d0
        gam(2)=0.d0
        gam(3)=0.d0
        rom(1)=-0.5d0
        rom(2)=0.d0
        rom(3)=0.d0
!m======================================================
      if(ismaster) then
        write(6,110) gam(1),rom(1)
  110   format(/,5x,'The time scheme is the Adams-Bashfort',4x, &
         'gam= ',f8.3,4x,'ro= ',f8.3)
      endif
     
!m======================================================                                 
      endif                                                             

      do ns=1,nsst
        alm(ns)=(gam(ns)+rom(ns))
      end do

      return
      end

