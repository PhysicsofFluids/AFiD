!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: solq3k.F90                                    !
!    CONTAINS: subroutine solq12k                         !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any the vertical direction, and updates it to    !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolveImpEqnUpdate_Z
      use param
      use local_arrays, only : q3,rhs
      use decomp_2d, only: xstart,xend
      implicit none
      real, dimension(n3) :: amkl,apkl,ackl, fkl
      real :: amkT(n3-1),apkT(n3-1)
      real :: appk(n3-2)
      real :: ackT(n3)
      integer :: jc,kc,info,ic
      integer :: ipkv(n3)
      real :: betadx,ackl_b

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
