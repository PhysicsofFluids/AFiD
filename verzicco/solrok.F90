!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: solrok.F90                                     !
!    CONTAINS: subroutine solrok                          !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for           !
!     temperature, and updates it to time t+dt            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine solrok
      use param
      use local_arrays, only : dens,rhs
      use decomp_2d, only: xstart,xend
      implicit none
      real, dimension(n3) :: amkl,apkl,ackl, fkl
      integer :: jc,kc,info,ipkv(n3),ic
      real :: betadx,ackl_b
      real :: amkT(n3-1),ackT(n3),apkT(n3-1),appk(n3-2)

!     Calculate the coefficients of the tridiagonal matrix
!     The coefficients are normalized to prevent floating
!     point errors.

      betadx=0.5d0*al*dt/pec

      amkl(1)=0.d0
      apkl(1)=0.d0
      ackl(1)=1.d0
      do kc=2,n3m
        ackl_b=1.0d0/(1.-ac3ssk(kc)*betadx)
        amkl(kc)=-am3ssk(kc)*betadx*ackl_b
        ackl(kc)=1.0d0
        apkl(kc)=-ap3ssk(kc)*betadx*ackl_b
      enddo
      amkl(n3)=0.d0
      apkl(n3)=0.d0
      ackl(n3)=1.d0

      amkT=amkl(2:n3)
      apkT=apkl(1:n3m)
      ackT=ackl(1:n3)

!     Call to LAPACK library to factor tridiagonal matrix.
!     No solving is done in this call.

      call dgttrf(n3,amkT,ackT,apkT,appk,ipkv,info)

      do ic=xstart(3),xend(3)
       do jc=xstart(2),xend(2)

!     Normalize RHS of equation

        fkl(1)= 0.d0
        do kc=2,n3m
         ackl_b=1.0/(1.-ac3ssk(kc)*betadx)
         fkl(kc)=rhs(kc,jc,ic)*ackl_b
        end do
        fkl(n3)= 0.d0
          
!     Solve equation using LAPACK library

        call dgttrs('N',n3,1,amkT,ackT,apkT,appk,ipkv,fkl,n3,info)
          
!      Update temperature field

        do kc=2,n3m
          dens(kc,jc,ic) = dens(kc,jc,ic) + fkl(kc)
        end do

       enddo
      end do

      return
      end
