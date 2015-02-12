!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVY.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVY             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the y (horizontal) dimension        !
!     and call the implicit solver                        !
!     After this routine, the velocity field in y has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateVY
      use param
      use local_arrays, only: vy,ruy,pr,rhs,dph
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,jmm,jc,ic
      integer :: kpp,kmm
      real    :: alre,udy
      real    :: amm,acc,app
      real    :: d33vy,dpx22


      alre=al/ren
      udy=dy*al

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vy,pr) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
!$OMP   SHARED(dy,al,ga,ro,alre,dt,dph) &
!$OMP   SHARED(udy,udx3m,rhs,ruy) &
!$OMP   PRIVATE(ic,jc,kc,kmm,kpp,jmm) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(d33vy,dpx22)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      jmm=jc-1
      do kc=1,nxm
      kmm=kmv(kc)
      kpp=kpv(kc)
      amm=am3sk(kc)
      acc=ac3sk(kc)
      app=ap3sk(kc)

!
!   33 second derivative of vy
!
            d33vy=vy(kpp,jc,ic)*app &
                 +vy(kc,jc,ic)*acc &
                 +vy(kmm,jc,ic)*amm

!
!   component of grad(pr) along 2 direction
!
            dpx22=(pr(kc,jc,ic)-pr(kc,jmm,ic))*udy

!
            rhs(kc,jc,ic)=(ga*dph(kc,jc,ic)+ro*ruy(kc,jc,ic) &
                          +alre*d33vy-dpx22)*dt

!m===========================================================
!
            ruy(kc,jc,ic)=dph(kc,jc,ic)
      enddo
      enddo
      enddo

      call SolveImpEqnUpdate_YZ(vy,rhs)
      
      return
      end
!
