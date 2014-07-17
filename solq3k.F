************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x3
c
      subroutine solq3k
      use param
      use local_arrays, only : q3,rhs
      use decomp_2d, only: xstart,xend
      implicit none
      real, dimension(m3) :: amkl,apkl,ackl, fkl
      real :: amkT(m3-1),apkT(m3-1)
      real :: appk(m3-2)
      real :: ackT(m3)
      integer :: jc,kc,info,ic
      integer :: ipkv(m3)
      real :: betadx,ackl_b

cm    dimension amkl(m3),apkl(m3),ackl(m3),fkl(m1,m3)
c  ********* compute the dq3* sweeping in the x3 direction
c
      betadx=beta*al

      amkl(1)=0.d0
      apkl(1)=0.d0
      ackl(1)=1.d0
      do kc=2,n3m
        ackl_b=1.0d0/(1.0d0-ac3ck(kc)*betadx)
        amkl(kc)=-am3ssk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3ssk(kc)*betadx*ackl_b
      enddo
      amkl(n3)=0.d0
      apkl(n3)=0.d0
      ackl(n3)=1.d0

      amkT=amkl(2:n3)
      apkT=apkl(1:(n3-1))
      ackT=ackl(1:n3)

      call dgttrf(n3,amkT,ackT,apkT,appk,ipkv,info)

      do ic=xstart(3),xend(3)
          do jc=xstart(2),xend(2)
            fkl(1)= 0.d0
          do kc=2,n3m
            ackl_b=1.0d0/(1.0d0-ac3ck(kc)*betadx)
            fkl(kc) = rhs(kc,jc,ic)*ackl_b
          enddo
            fkl(n3)= 0.d0
          
          call dgttrs('N',n3,1,amkT,ackT,apkT,appk,ipkv,fkl,n3,info)

          do kc=2,n3m
            q3(kc,jc,ic)=q3(kc,jc,ic) + fkl(kc)
          enddo
          enddo
      end do

      return
      end
