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
      use local_arrays, only: q2,q3,temp,q1
      use decomp_2d, only: xstart,xend
      use mpih
      implicit none
      integer :: j,k,i
      real :: xxx,yyy,eps

      eps=1.0d0
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=1,nxm
           yyy=zm(k) 
           q1(k,j,i)=0.0d0
           yyy=zm(k) 
           xxx=rc(j)            
           q2(k,j,i)=(2.0d0*yyy-6.0d0*yyy**2+4.0d0*yyy**3) &
     &                  *sin(3*xxx)*eps

           yyy=zz(k)          
           xxx=rm(j)
           q3(k,j,i)=-yyy**2*(1.0d0-yyy)**2*cos(3.1*xxx)*eps

         enddo
        enddo
      enddo

      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=2,nxm
             temp(k,j,i)= tempbp(j,i) - (tempbp(j,i) - temptp(j,i)) &
                         *zz(k)
           enddo
          end do 
        end do
      

        temp(1,:,:)=1.0d0
        temp(nx,:,:)=0.0d0

      return                                                            
      end                                                               
