!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: invtr2.F90                                     !
!    CONTAINS: subroutine invtr2                          !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the second horizontal dimension     !
!     and call the implicit solver                        !
!     After this routine, the velocity field in x2 has    !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateVY
      use param
      use local_arrays, only: q2,ru2,pr,rhs,dph
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,jmm,jc,ic
      integer :: kpp,kmm
      real    :: alre,udx2
      real    :: amm,acc,app
      real    :: d33q2,dpx22


      alre=al/ren
      udx2=dx2*al

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,q2,pr) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
!$OMP   SHARED(dx2,al,ga,ro,alre,dt,dph) &
!$OMP   SHARED(udx2,udx3m,rhs,ru2) &
!$OMP   PRIVATE(ic,jc,kc,kmm,kpp,jmm) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(d33q2,dpx22)
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
!   33 second derivative of q2
!
            d33q2=q2(kpp,jc,ic)*app &
                 +q2(kc,jc,ic)*acc &
                 +q2(kmm,jc,ic)*amm

!
!   component of grad(pr) along 2 direction
!
            dpx22=(pr(kc,jc,ic)-pr(kc,jmm,ic))*udx2

!
            rhs(kc,jc,ic)=(ga*dph(kc,jc,ic)+ro*ru2(kc,jc,ic) &
                          +alre*d33q2-dpx22)*dt

!m===========================================================
!
            ru2(kc,jc,ic)=dph(kc,jc,ic)
      enddo
      enddo
      enddo

      call SolveImpEqnUpdate_XY(q2,rhs)
      
      return
      end
!
