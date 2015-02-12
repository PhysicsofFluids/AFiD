!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVX.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVX             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the X (vertical) direction and call !
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
      real    :: amm,acc,app,dxp,dxxvx

      alre=al/ren

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vx,pr) &
!$OMP   SHARED(kmv,kpv,am3ck,ac3ck,ap3ck) &
!$OMP   SHARED(al,ga,ro,alre,dt,qcap) &
!$OMP   SHARED(udx3c,rhs,rux) &
!$OMP   PRIVATE(ic,jc,kc,km,kp) &
!$OMP   PRIVATE(amm,acc,app,udx3) &
!$OMP   PRIVATE(dxxvx,dxp)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      do kc=2,nxm
      km=kc-1
      kp=kc+1
      udx3 = al*udx3c(kc)
      amm=am3ck(kc)
      acc=ac3ck(kc)
      app=ap3ck(kc)

!   Second derivative in x-direction of vx
!
            dxxvx=vx(kp,jc,ic)*app &
                +vx(kc,jc,ic)*acc &
                +vx(km,jc,ic)*amm

!  component of grad(pr) along x direction
!
            dxp=(pr(kc,jc,ic)-pr(km,jc,ic))*udx3

!    Calculate right hand side of Eq. 5 (VO96)
!
            rhs(kc,jc,ic)=(ga*qcap(kc,jc,ic)+ro*rux(kc,jc,ic) &
                          +alre*dxxvx-dxp)*dt 

!    Store the non-linear terms for the calculation of 
!    the next timestep

            rux(kc,jc,ic)=qcap(kc,jc,ic)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

!  Solve equation and update velocity

      call SolveImpEqnUpdate_X

!  Set boundary conditions on the vertical velocity at top
!  and bottom plates. This seems necessary.

      vx(1,:,:)=0.0d0
      vx(nx,:,:)=0.0d0

      return
      end
!
