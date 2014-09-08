!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: hdnl1.F90                                      !
!    CONTAINS: subroutine hdnl1                           !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the first horizontal dimension.     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hdnl1
      use param
      use local_arrays, only: q2,q3,dph,q1,dq
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,kp,jpp,jmm,jc,ic,imm,ipp
      integer :: kmm,kpp
      real    :: h11,h12,h13,udx1,udx2
      real    :: udx1q,udx2q
      real    :: d11q1,d22q1

!
      udx1q=dx1q/ren
      udx2q=dx2q/ren

      udx1=dx1*0.25
      udx2=dx2*0.25

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,q1,q2,q3,dx1,dx2,udx3m) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udx1) &
!$OMP   SHARED(udx2,udx1q,udx2q,dq,n3m) &
!$OMP   PRIVATE(ic,jc,kc,imm,ipp,kmm,kp,kpp) &
!$OMP   PRIVATE(jmm,jpp) &
!$OMP   PRIVATE(h11,h12,h13,d11q1,d22q1)

      do ic=xstart(3),xend(3)
       imm=ic-1
       ipp=ic+1
       do jc=xstart(2),xend(2)
        jmm=jc-1
        jpp=jc+1
        do kc=1,n3m
         kmm=kmv(kc)
         kpp=kpv(kc)
         kp=kc+1
      
!     q1 q1 term
!
!
!                 d  q_t q_t 
!                ------------
!                 d   t      
!
      h11=( (q1(kc,jc,ipp)+q1(kc,jc,ic)) &
           *(q1(kc,jc,ipp)+q1(kc,jc,ic)) &
           -(q1(kc,jc,imm)+q1(kc,jc,ic)) &
           *(q1(kc,jc,imm)+q1(kc,jc,ic)) &
          )*udx1

!     q1 q2 term
!
!
!                 d  q_t q_r 
!                ------------
!                 d   r      
!
      h12=( (q2(kc,jpp,ic)+q2(kc,jpp,imm)) &
           *(q1(kc,jpp,ic)+q1(kc,jc,ic)) &
           -(q2(kc,jc,ic)+q2(kc,jc,imm)) &
           *(q1(kc,jc,ic)+q1(kc,jmm,ic)) &
          )*udx2
!
!     q1 q3 term
!
!
!                 d  q_t q_x 
!                -----------
!                 d   x      
!
      h13=((q3(kp,jc,ic)+q3(kp,jc,imm))*(q1(kpp,jc,ic)+q1(kc,jc,ic)) &
          -(q3(kc,jc,ic)+q3(kc,jc,imm))*(q1(kc,jc,ic)+q1(kmm,jc,ic)) &
          )*udx3m(kc)*0.25d0
!
!
!
!   11 second derivative of q1
!
            d11q1=(q1(kc,jc,ipp) &
                  -2.0*q1(kc,jc,ic) &
                  +q1(kc,jc,imm))*udx1q
!
!   22 second derivative of q1
!
            d22q1=(q1(kc,jpp,ic) &
                  -2.0*q1(kc,jc,ic) &
                  +q1(kc,jmm,ic))*udx2q

!
        dq(kc,jc,ic)=-(h11+h12+h13)+d22q1+d11q1
!
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
!
