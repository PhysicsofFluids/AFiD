!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ImplicitAndUpdateVZ.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVZ             !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the z (horizontal) dimension        !
!     and call the implicit solver.                       !
!     After this routine, the velocity field in z has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateVZ
      use param
      use local_arrays, only: vz,dq,ruz,rhs,pr
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,jc,ic,imm
      integer :: kmm,kpp
      real    :: alre,amm,acc,app
      real    :: d33vz,dpx11,udz

      alre=al/ren
      udz=dz*al

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vz,pr) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
!$OMP   SHARED(dz,al,ga,ro,alre,dt,dq) &
!$OMP   SHARED(udx3m,rhs,ruz) &
!$OMP   PRIVATE(ic,jc,kc,imm,kmm,kpp) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(d33vz,dpx11)
      do ic=xstart(3),xend(3)
      imm=ic-1
      do jc=xstart(2),xend(2)
      do kc=1,nxm
      kmm=kmv(kc)
      kpp=kpv(kc)
      amm=am3sk(kc)
      acc=ac3sk(kc)
      app=ap3sk(kc)

!
!   33 second derivative of vz
!
            d33vz=vz(kpp,jc,ic)*app &
                 +vz(kc,jc,ic)*acc &
                 +vz(kmm,jc,ic)*amm
      
!
!   component of grad(pr) along 2 direction
!
            dpx11=(pr(kc,jc,ic)-pr(kc,jc,imm))*dz*al
!
!
!
            rhs(kc,jc,ic)=(ga*dq(kc,jc,ic)+ro*ruz(kc,jc,ic) &
                          +alre*d33vz-dpx11)*dt

       

!m===========================================================
!
            ruz(kc,jc,ic)=dq(kc,jc,ic)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO


      call SolveImpEqnUpdate_YZ(vz,rhs)


      return
      end
!
