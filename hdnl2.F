c***********************************************************************
      subroutine hdnl2
      use param
      use local_arrays, only: q2,q3,q1,dph
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,kp,jpp,jmm,jc,ic,imm,ipp
      integer :: kpp,kmm
      real    :: udx1q,udx2q
      real    :: h22,h23,udx1,udx2,h21
      real    :: d11q2,d22q2


      udx1q=dx1q/ren
      udx2q=dx2q/ren
      udx1=dx1*0.25
      udx2=dx2*0.25
c
c     h term for the q2 momentum equation at i+1/2,j,k+1/2
c
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(xstart,xend,q1,q2,q3,dx1,dx2)
!$OMP$ SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udx1)
!$OMP$ SHARED(udx2,udx1q,udx2q,udx3m,dph,n3m)
!$OMP$ PRIVATE(ic,jc,kc,imm,ipp,kmm,kp,kpp)
!$OMP$ PRIVATE(jmm,jpp)
!$OMP$ PRIVATE(h21,h22,h23,d11q2,d22q2)
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

c     q1 q2 term
c
c
c                 d  q_t q_r 
c                ------------
c                 d   t      
c
      h21=( (q2(kc,jc,ipp)+q2(kc,jc,ic))
     %     *(q1(kc,jc,ipp)+q1(kc,jmm,ipp))
     %     -(q2(kc,jc,ic)+q2(kc,jc,imm))
     %     *(q1(kc,jc,ic)+q1(kc,jmm,ic))
     %    )*udx1
      
c     q2 q2 term
c
c
c                 d  q_r q_r 
c                ------------
c                 d   r      
c
      h22=( (q2(kc,jpp,ic)+q2(kc,jc,ic))
     %     *(q2(kc,jpp,ic)+q2(kc,jc,ic))
     %     -(q2(kc,jmm,ic)+q2(kc,jc,ic))
     %     *(q2(kc,jmm,ic)+q2(kc,jc,ic))
     %    )*udx2
c
c     q2 q3 term
c
c
c                 d  q_x q_r 
c                -----------
c                 d   x      
c
      h23=((q3(kp,jc,ic)+q3(kp,jmm,ic))*(q2(kpp,jc,ic)+q2(kc,jc,ic))
     %    -(q3(kc,jc,ic)+q3(kc,jmm,ic))*(q2(kc,jc,ic)+q2(kmm,jc,ic))
     %    )*udx3m(kc)*0.25d0
c
c
c
c   11 second derivative of q2
c
            d11q2=(q2(kc,jc,ipp)
     %            -2.0*q2(kc,jc,ic)
     %            +q2(kc,jc,imm))*udx1q

c
c   22 second derivative of q2
c
            d22q2=(q2(kc,jpp,ic)
     %            -2.0*q2(kc,jc,ic)
     %            +q2(kc,jmm,ic))*udx2q

            dph(kc,jc,ic)=-(h21+h22+h23)+d11q2+d22q2
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
c
