!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: invtr1.F90                                     !
!    CONTAINS: subroutine invtr1                          !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the first horizontal dimension      !
!     and call the implicit solver.                       !
!     After this routine, the velocity field in x1 has    !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateVX
      use param
      use local_arrays, only: q1,dq,ru1,rhs,pr
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,jc,ic,imm
      integer :: kmm,kpp
      real    :: alre,amm,acc,app
      real    :: d33q1,dpx11

!
      alre=al/ren


!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,n3m,q1,pr) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk) &
!$OMP   SHARED(dx1,al,ga,ro,alre,dt,dq) &
!$OMP   SHARED(udx3m,rhs,ru1) &
!$OMP   PRIVATE(ic,jc,kc,imm,kmm,kpp) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(d33q1,dpx11)
      do ic=xstart(3),xend(3)
      imm=ic-1
      do jc=xstart(2),xend(2)
      do kc=1,n3m
      kmm=kmv(kc)
      kpp=kpv(kc)
      amm=am3sk(kc)
      acc=ac3sk(kc)
      app=ap3sk(kc)

!
!   33 second derivative of q1
!
            d33q1=q1(kpp,jc,ic)*app &
                 +q1(kc,jc,ic)*acc &
                 +q1(kmm,jc,ic)*amm
      
!
!   component of grad(pr) along 2 direction
!
            dpx11=(pr(kc,jc,ic)-pr(kc,jc,imm))*dx1*al
!
!
!
            rhs(kc,jc,ic)=(ga*dq(kc,jc,ic)+ro*ru1(kc,jc,ic) &
                          +alre*d33q1-dpx11)*dt

       

!m===========================================================
!
            ru1(kc,jc,ic)=dq(kc,jc,ic)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO


      call SolveImpEqnUpdate_XY(q1,rhs)


      return
      end
!
