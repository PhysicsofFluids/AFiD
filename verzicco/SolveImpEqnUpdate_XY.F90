!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: solq12k.F90                                    !
!    CONTAINS: subroutine solq12k                         !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for velocity  !
!     in any of the horizontal directions, and updates    !
!     it to time t+dt                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine SolveImpEqnUpdate_XY(q,rhs)
      use param
      use decomp_2d, only: xstart,xend
      implicit none
      real,intent(inout) :: q(1:n3,xstart(2)-1:xend(2)+1, &
     &      xstart(3)-1:xend(3)+1)
      real,intent(inout) :: rhs(1:n3,xstart(2):xend(2), &
     &      xstart(3):xend(3))
      real, dimension(n3) :: amkl,apkl,ackl,fkl
      integer :: jc,kc,info,ipkv(n3m),ic
      real :: betadx,ackl_b
      real :: amkT(n3m-1),ackT(n3m),apkT(n3m-1),appk(n3-3)

      betadx=beta*al

      do kc=1,n3m
        ackl_b=1.0d0/(1.0d0-ac3sk(kc)*betadx)
        amkl(kc)=-am3sk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3sk(kc)*betadx*ackl_b
      enddo

      amkT=amkl(2:n3m)
      apkT=apkl(1:(n3m-1))
      ackT=ackl(1:n3m)

      call dgttrf(n3m,amkT,ackT,apkT,appk,ipkv,info)

      do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
          do kc=1,n3m
            ackl_b=1.0d0/(1.-ac3sk(kc)*betadx)
            fkl(kc)=rhs(kc,jc,ic)*ackl_b
          end do

          call dgttrs('N',n3m,1,amkT,ackT,apkT,appk,ipkv,fkl,n3m,info)

          do kc=1,n3m
            q(kc,jc,ic)=q(kc,jc,ic) + fkl(kc)
          end do
        end do
      end do

      return
      end
