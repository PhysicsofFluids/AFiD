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
      real, dimension(nx) :: amkl,apkl,ackl,fkl
      integer :: jc,kc,info,ipkv(nxm),ic
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

      do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
          do kc=1,nxm
            ackl_b=1.0d0/(1.-ac3sk(kc)*betadx)
            fkl(kc)=rhs(kc,jc,ic)*ackl_b
          end do

          call dgttrs('N',nxm,1,amkT,ackT,apkT,appk,ipkv,fkl,nxm,info)

          do kc=1,nxm
            q(kc,jc,ic)=q(kc,jc,ic) + fkl(kc)
          end do
        end do
      end do

      return
      end
