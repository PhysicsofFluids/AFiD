!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: hdnl3.F90                                      !
!    CONTAINS: subroutine hdnl3                           !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the vertical dimension (x3)         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hdnl3
      use param
      use local_arrays, only: q1,q2,q3,dens,qcap
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,kc
      integer :: km,kp,jmm,jpp,ic,imm,ipp
      real    :: h32,h33,h31
      real    :: udx1,udx2,densit
      real    :: udx1q,udx2q
      real    :: dq31,dq32

      udx1=dx1*0.25
      udx2=dx2*0.25
      udx1q=dx1q/ren
      udx2q=dx2q/ren

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,n3m,q1,q2,q3,dx1,dx2) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udx1) &
!$OMP   SHARED(udx2,udx1q,udx2q,udx3c,qcap,dens) &
!$OMP   PRIVATE(ic,jc,kc,imm,ipp,km,kp) &
!$OMP   PRIVATE(jmm,jpp,densit) &
!$OMP   PRIVATE(h31,h32,h33,dq31,dq32)
      do ic=xstart(3),xend(3)
      imm=ic-1
      ipp=ic+1
      do jc=xstart(2),xend(2)
      jmm=jc-1
      jpp=jc+1
      do kc=2,n3m
      km=kc-1
      kp=kc+1
!
!    q3 q1 term
!
!
!                d  q_x q_t 
!             -----------
!                d   t      
!
!
      h31=(((q1(kc,jc,ipp)+q1(km,jc,ipp)) &
           *(q3(kc,jc,ipp)+q3(kc,jc,ic))) &
          -((q1(kc,jc,ic)+q1(km,jc,ic)) &
           *(q3(kc,jc,ic)+q3(kc,jc,imm))))*udx1
!
!    q3 q2 term
!
!
!                d  q_x q_r 
!             -----------
!                d   r      
!
      h32=(((q2(kc,jpp,ic)+q2(km,jpp,ic)) &
           *(q3(kc,jpp,ic)+q3(kc,jc,ic))) &
          -((q2(kc,jc,ic)+q2(km,jc,ic)) &
           *(q3(kc,jc,ic)+q3(kc,jmm,ic))))*udx2
!
!    q3 q3 term
!
!
!                 d  q_x q_x 
!                -----------
!                 d   x      
!
      h33=((q3(kp,jc,ic)+q3(kc,jc,ic))*(q3(kp,jc,ic)+q3(kc,jc,ic)) &
          -(q3(kc,jc,ic)+q3(km,jc,ic))*(q3(kc,jc,ic)+q3(km,jc,ic)) &
          )*udx3c(kc)*0.25d0
!
!  add the buoyancy term
!
          densit=dens(kc,jc,ic)

!
!   11 second derivatives of q3
!
            dq31=(q3(kc,jc,imm) &
                 -2.0*q3(kc,jc,ic) &
                 +q3(kc,jc,ipp))*udx1q
!
!   22 second derivatives of q3
!
            dq32=(q3(kc,jmm,ic) &
                 -2.0*q3(kc,jc,ic) &
                 +q3(kc,jpp,ic))*udx2q


          qcap(kc,jc,ic) =-(h31+h32+h33)+dq31+dq32+densit
            
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
!
