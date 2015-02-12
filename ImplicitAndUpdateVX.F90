!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVX.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVX             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the vertical direction and call     !
!     the implicit solver. After this routine, the        !
!     vertical velocity has been updated to the new       !
!     timestep                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateVX
      use param
      use local_arrays, only: vx,rhs,rux,qcap,pr
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,kc
      integer :: km,kp,ic
      real    :: alre,udx3
      real    :: amm,acc,app,dpx33,dq33

      alre=al/ren

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vx,pr) &
!$OMP   SHARED(kmv,kpv,am3ck,ac3ck,ap3ck) &
!$OMP   SHARED(al,ga,ro,alre,dt,qcap) &
!$OMP   SHARED(udx3c,rhs,rux) &
!$OMP   PRIVATE(ic,jc,kc,km,kp) &
!$OMP   PRIVATE(amm,acc,app,udx3) &
!$OMP   PRIVATE(dq33,dpx33)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      do kc=2,nxm
      km=kc-1
      kp=kc+1
      udx3 = al*udx3c(kc)
      amm=am3ck(kc)
      acc=ac3ck(kc)
      app=ap3ck(kc)

!   33 second derivatives of vx
!
            dq33=vx(kp,jc,ic)*app &
                +vx(kc,jc,ic)*acc &
                +vx(km,jc,ic)*amm

!  component of grad(pr) along x3 direction
!
            dpx33=(pr(kc,jc,ic)-pr(km,jc,ic))*udx3
!m=======================================================     
            rhs(kc,jc,ic)=(ga*qcap(kc,jc,ic)+ro*rux(kc,jc,ic) &
                          +alre*dq33-dpx33)*dt 
!m=======================================================
!
!  updating of the non-linear terms
!
            rux(kc,jc,ic)=qcap(kc,jc,ic)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO


      call SolveImpEqnUpdate_X

      vx(1,:,:)=0.0d0
      vx(nx,:,:)=0.0d0

      return
      end
!
