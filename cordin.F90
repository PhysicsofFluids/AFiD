!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: cordin.F90                                     !
!    CONTAINS: subroutine cordin                          !
!                                                         ! 
!    PURPOSE: Initialization routine. Calcuates the       !
!     coordinates of points in the three dimensional grid !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine cordin
      use param                                                 
      use mpih
      use decomp_2d, only: nrank
      implicit none
      integer  :: j,k,n3mo,nclip,i
      real, dimension(1:m3) :: etaz
      real, dimension(1:m3+500) :: etazm
      real :: tstr3, z2dp
      real :: x1,x2,x3,etain,delet

!
!     RADIAL COORDINATE DEFINITION
!
      rint = 0.d0
!
!     UNIFORM GRID
!
        do  i=1,n1
         x1=real(i-1)/real(n1m)
         tc(i)= rext*x1
        end do
      do i=1,n1m
        tm(i)=(tc(i)+tc(i+1))*0.5d0
      end do

        do  j=1,n2
         x2=real(j-1)/real(n2m)
         rc(j)= rext2*x2
        end do
      do j=1,n2m
        rm(j)=(rc(j)+rc(j+1))*0.5d0
      end do
!
!     AXIAL COORDINATE DEFINITION
!
!
!     UNIFORM GRID
!
!       write(6,*) istr3,alx3

      tstr3=tanh(str3)

      if (istr3.eq.0) then
        do k=1,n3
          x3=real(k-1)/real(n3m)
          etaz(k)=alx3*x3
          zz(k)=etaz(k)
        enddo
      endif


!       
!      CLUSTERING AT THE EXTERNAL RADIAL WALL 
!                       and  
!             CLUSTERING AT THE AXIS 
!      

        if (istr3.eq.4) then
         zz(1)=0.0d0
         do k=2,n3
          z2dp=float(2*k-n3-1)/float(n3m)
          zz(k)=(1+tanh(str3*z2dp)/tstr3)*0.5*alx3
          if(zz(k).lt.0.or.zz(k).gt.alx3)then
           write(*,*)'Forza la griglia: ','zc(',k,')=',zz(k)
           stop
          endif
         end do
        end if



      if(istr3.eq.6) then
      nclip = int(str3)
      n3mo = n3+nclip+nclip
      do k=1,n3mo
        etazm(k)=+cos(pi*(float(k)-0.5)/float(n3mo))
      end do
      do k=1,n3
        etaz(k)=etazm(k+nclip)
      end do
      delet = etaz(1)-etaz(n3)
      etain = etaz(1)
      do k=1,n3
        etaz(k)=etaz(k)/(0.5*delet)
      end do
      zz(1) = 0.
      do k=2,n3m
        zz(k) = alx3*(1.-etaz(k))*0.5
      end do
      zz(n3) = alx3
      endif
      
!m-----------------------------------------
!
!     STAGGERED COORDINATES AND
!     METRIC QUANTITIES
!
      do k=1,n3m
        zm(k)=(zz(k)+zz(k+1))*0.5d0
        g3rm(k)=(zz(k+1)-zz(k))*dx3
      enddo
      do k=2,n3m
        g3rc(k)=(zz(k+1)-zz(k-1))*dx3*0.5d0
      enddo
      g3rc(1)=(zz(2)-zz(1))*dx3
      g3rc(n3)= (zz(n3)-zz(n3m))*dx3
!
!     WRITE GRID INFORMATION
!
      do k=1,n3m
        udx3m(k) = dx3/g3rm(k)
        udx3c(k) = dx3/g3rc(k)
      end do
      udx3c(n3) = dx3/g3rc(n3)
!m====================================================
      if(nrank.eq.0) then
      open(unit=78,file='axicor.out',status='unknown')
      do k=1,n3
        write(78,345) k,zz(k),zm(k),g3rc(k),g3rm(k)
      end do
      close(78)
 345  format(i4,4(2x,e23.15))
!m===================================================
!
!     QUANTITIES FOR DERIVATIVES
!
      open(unit=78,file='fact3.out',status='unknown')
      do k=1,n3m
        write(78,*) k,udx3m(k),udx3c(k)
      end do
        write(78,*) n3,udx3m(n3m),udx3c(n3)
      close(78)

      endif

      return                                                            
      end                                                               
