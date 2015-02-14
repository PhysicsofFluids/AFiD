!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CreateInitialConditions.F90                    !
!    CONTAINS: subroutine CreateInitialConditions         !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets initial        !
!     conditions for velocity and temperature             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CreateInitialConditions
      use param
      use local_arrays, only: vy,vx,temp,vz
      use decomp_2d, only: xstart,xend
      use mpih
      implicit none
      integer :: j,k,i
      real :: xxx,yyy,eps

      eps=1.0d0
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=1,nxm
           vz(k,j,i)=0.0d0
           yyy=xm(k) 
           xxx=yc(j)            
           vy(k,j,i)=(2.0d0*yyy-6.0d0*yyy**2+4.0d0*yyy**3) &
     &                  *sin(3*xxx)*eps

           yyy=xc(k)          
           xxx=ym(j)
           vx(k,j,i)=-yyy**2*(1.0d0-yyy)**2*cos(3.1*xxx)*eps

         enddo
        enddo
      enddo

      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=2,nxm
             temp(k,j,i)= tempbp(j,i) - (tempbp(j,i) - temptp(j,i)) &
                         *xc(k)
           enddo
          end do 
        end do
      

        temp(1,:,:)=1.0d0
        temp(nx,:,:)=0.0d0

      return                                                            
      end                                                               
