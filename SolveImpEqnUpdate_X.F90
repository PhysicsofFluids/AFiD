!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_X.F90                        !
!    CONTAINS: subroutine SolveImpEqnUpdate_X             !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any the vertical direction, and updates it to    !
!     time t+dt                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolveImpEqnUpdate_X
      use param
      use local_arrays, only : vx,rhs
      use decomp_2d, only: xstart,xend
      implicit none
      real, dimension(nx) :: amkl,apkl,ackl
      real :: amkT(nx-1),apkT(nx-1)
      real :: appk(nx-2)
      real :: ackT(nx)
      integer :: jc,kc,info,ic,nrhs
      integer :: ipkv(nx)
      real :: betadx,ackl_b

      betadx=beta*al

      amkl(1)=0.d0
      apkl(1)=0.d0
      ackl(1)=1.d0
      do kc=2,nxm
        ackl_b=1.0d0/(1.0d0-ac3ck(kc)*betadx)
        amkl(kc)=-am3ssk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3ssk(kc)*betadx*ackl_b
      enddo
      amkl(nx)=0.d0
      apkl(nx)=0.d0
      ackl(nx)=1.d0

      amkT=amkl(2:nx)
      apkT=apkl(1:(nx-1))
      ackT=ackl(1:nx)

      call dgttrf(nx,amkT,ackT,apkT,appk,ipkv,info)

      nrhs=(xend(3)-xstart(3)+1)*(xend(2)-xstart(2)+1)
      do ic=xstart(3),xend(3)
         do jc=xstart(2),xend(2)
            do kc=2,nxm
              ackl_b=1.0/(1.0-ac3ck(kc)*betadx)
              rhs(kc,jc,ic)=rhs(kc,jc,ic)*ackl_b
            end do
         end do
      end do

      call dgttrs('N',nx,nrhs,amkT,ackT,apkT,appk,ipkv,rhs,nx,info)

       do ic=xstart(3),xend(3)
         do jc=xstart(2),xend(2)
            do kc=2,nxm
              vx(kc,jc,ic)=vx(kc,jc,ic) + rhs(kc,jc,ic)
             end do
          end do
      end do

      return
      end
