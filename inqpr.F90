!***********************************************************************
!                                                                      *
!                       INITIAL CONDITION                              *
!                                                                      *
!***********************************************************************
      subroutine inqpr
      use param
      use local_arrays, only: q2,q3,dens,q1
      use decomp_2d, only: xstart,xend
      use mpih
      implicit none
      integer :: j,k,i
      real :: xxx,yyy,eps

      eps=1.0d0
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=1,n3m
           yyy=zm(k) 
!           q1(k,j,i)=0.1*yyy*(1-yyy)
           q1(k,j,i)=0.0d0
           yyy=zm(k) 
           xxx=rc(j)            
           q2(k,j,i)=(2.0d0*yyy-6.0d0*yyy**2+4.0d0*yyy**3) &
     &                  *sin(3*xxx)*eps

           yyy=zz(k)          
           xxx=rm(j)
           q3(k,j,i)=-yyy**2*(1.0d0-yyy)**2*cos(3.1*xxx)*eps

!        q2(k,j,i)=0.0d0
!        q3(k,j,i)=0.0d0
         enddo
        enddo
      enddo

      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=2,n3m
             dens(k,j,i)= denbs(j,i) - (denbs(j,i) - denbn(j,i)) &
                         *zz(k)
!             dens(k,j,i)= (i+j)*k
           enddo
          end do 
        end do
      

        dens(1,:,:)=1.0d0
        dens(n3,:,:)=0.0d0

#ifdef SERIAL_DEBUG        
        mck1=0.0d0
        mck2=0.0d0
        mck3=0.0d0
        mck4=0.0d0
        do k=1,n3m
        do j=1,n2m
        do i=1,n1m
        mck1=mck1+q1(k,j,i)
        mck2=mck2+q2(k,j,i)
        mck3=mck3+q3(k,j,i)
        mck4=mck3+dens(k,j,i)
        enddo
        enddo
        enddo
        write(*,*) 'q1cksum= ',mck1
        write(*,*) 'q2cksum= ',mck2
        write(*,*) 'q3cksum= ',mck3
        write(*,*) 'denscksum= ',mck4
        write(*,*) 'starting tsch'
#endif
      return                                                            
      end                                                               
!
