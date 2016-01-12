!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdate_YZ.F90                       !
!    CONTAINS: subroutine SolveImpEqnUpdate_YZ            !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any of the horizontal directions, and updates    !
!     it to time t+dt                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine SolveImpEqnUpdate_YZ(q,rhs)
      use param
      use decomp_2d, only: xstart,xend
      implicit none
      real,intent(inout) :: q(1:nx,xstart(2)-lvlhalo:xend(2)+lvlhalo, &
     &      xstart(3)-lvlhalo:xend(3)+lvlhalo)
      real,intent(inout) :: rhs(1:nx,xstart(2):xend(2), &
     &      xstart(3):xend(3))
      real, dimension(nx) :: amkl,apkl,ackl
      integer :: jc,kc,info,ipkv(nxm),ic,nrhs
      real :: betadx,ackl_b
      real :: amkT(nxm-1),ackT(nxm),apkT(nxm-1),appk(nx-3)

      betadx=beta*al

      do kc=1,nxm
        ackl_b=1.0d0/(1.0d0-ac3sk(kc)*betadx)
        amkl(kc)=-am3sk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3sk(kc)*betadx*ackl_b
      enddo

      amkT=amkl(2:nxm)
      apkT=apkl(1:(nxm-1))
      ackT=ackl(1:nxm)

      call dgttrf(nxm,amkT,ackT,apkT,appk,ipkv,info)

      nrhs=(xend(3)-xstart(3)+1)*(xend(2)-xstart(2)+1)
       do ic=xstart(3),xend(3)
         do jc=xstart(2),xend(2)
            do kc=1,nxm
               ackl_b=1.0/(1.0-ac3sk(kc)*betadx)
               rhs(kc,jc,ic)=rhs(kc,jc,ic)*ackl_b
             end do
          end do
      end do

      call dgttrs('N',nxm,nrhs,amkT,ackT,apkT,appk,ipkv,rhs,nx,info)

       do ic=xstart(3),xend(3)
         do jc=xstart(2),xend(2)
            do kc=1,nxm
              q(kc,jc,ic)=q(kc,jc,ic) + rhs(kc,jc,ic)
             end do
          end do
      end do

      return
      end
