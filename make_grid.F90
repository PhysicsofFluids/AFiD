!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: make_grid.F90                                  !
!    CONTAINS: subroutine make_grid                       !
!                                                         ! 
!    PURPOSE: Compute the indices, grid, grid metrics     !
!     and coefficients for differentiation                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine make_grid
      use param
      use decomp_2d, only: nrank
      implicit none

      real :: x1,x2,x3
      real :: a33, a33m, a33p
      real :: delet, etain, tstr3
      real :: z2dp

      real, dimension(m3+500) :: etaz, etazm

      integer :: i, j, kc, km, kp
      integer :: n3mo, nclip

      do kc=1,n3m
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=kc
        if(kc.eq.n3m) kpv(kc)=kc
      end do

      do kc=1,n3m
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
        kup(kc)=1-kpc(kc)
        kum(kc)=1-kmc(kc)
      end do


!
!     UNIFORM (HORIZONTAL DIRECTIONS) GRID
!
       do  i=1,n1
        x1=real(i-1)/real(n1m)
        tc(i)= rext*x1
       end do

       do i=1,n1m
         tm(i)=(tc(i)+tc(i+1))*0.5d0
       end do

       do j=1,n2
        x2=real(j-1)/real(n2m)
        rc(j)= rext2*x2
       end do

       do j=1,n2m
        rm(j)=(rc(j)+rc(j+1))*0.5d0
       end do

!
!     VERTICAL COORDINATE DEFINITION
!
!     OPTION 0: UNIFORM CLUSTERING
!

      if (istr3.eq.0) then
        do kc=1,n3
          x3=real(kc-1)/real(n3m)
          etaz(kc)=alx3*x3
          zz(kc)=etaz(kc)
        enddo
      endif

!
!     OPTION 4: HYPERBOLIC TANGENT-TYPE CLUSTERING
!

        tstr3=tanh(str3)

        if (istr3.eq.4) then
         zz(1)=0.0d0
         do kc=2,n3
          z2dp=float(2*kc-n3-1)/float(n3m)
          zz(kc)=(1+tanh(str3*z2dp)/tstr3)*0.5*alx3
          if(zz(kc).lt.0.or.zz(kc).gt.alx3)then
           write(*,*)'Forza la griglia: ','zc(',kc,')=',zz(kc)
           stop
          endif
         end do
        end if

!
!     OPTION 6: CLIPPED CHEBYCHEV-TYPE CLUSTERING
!


      if(istr3.eq.6) then
      nclip = int(str3)
      n3mo = n3+nclip+nclip
      do kc=1,n3mo
        etazm(kc)=+cos(pi*(float(kc)-0.5)/float(n3mo))
      end do
      do kc=1,n3
        etaz(kc)=etazm(kc+nclip)
      end do
      delet = etaz(1)-etaz(n3)
      etain = etaz(1)
      do kc=1,n3
        etaz(kc)=etaz(kc)/(0.5*delet)
      end do
      zz(1) = 0.
      do kc=2,n3m
        zz(kc) = alx3*(1.-etaz(kc))*0.5
      end do
      zz(n3) = alx3
      endif
      
!m-----------------------------------------
!
!     METRIC FOR UNIFORM DIRECTIONS
!

      dx1=real(n1m)/rext
      dx2=real(n2m)/rext2
      dx3=real(n3m)/alx3

      dx1q=dx1*dx1                                                      
      dx2q=dx2*dx2                                                      
      dx3q=dx3*dx3                                                      

!
!     STAGGERED COORDINATES AND
!     METRIC QUANTITIES FOR NON-UNIFORM 
!     DIRECTIONS
!

      do kc=1,n3m
        zm(kc)=(zz(kc)+zz(kc+1))*0.5d0
        g3rm(kc)=(zz(kc+1)-zz(kc))*dx3
      enddo
      do kc=2,n3m
        g3rc(kc)=(zz(kc+1)-zz(kc-1))*dx3*0.5d0
      enddo
      g3rc(1)=(zz(2)-zz(1))*dx3
      g3rc(n3)= (zz(n3)-zz(n3m))*dx3
!
!     WRITE GRID INFORMATION
!
      do kc=1,n3m
        udx3m(kc) = dx3/g3rm(kc)
        udx3c(kc) = dx3/g3rc(kc)
      end do
      udx3c(n3) = dx3/g3rc(n3)
!m====================================================
      if(nrank.eq.0) then
      open(unit=78,file='axicor.out',status='unknown')
      do kc=1,n3
        write(78,345) kc,zz(kc),zm(kc),g3rc(kc),g3rm(kc)
      end do
      close(78)
 345  format(i4,4(2x,e23.15))
!m===================================================
!
!     QUANTITIES FOR DERIVATIVES
!
      open(unit=78,file='fact3.out',status='unknown')
      do kc=1,n3m
        write(78,*) kc,udx3m(kc),udx3c(kc)
      end do
        write(78,*) n3,udx3m(n3m),udx3c(n3)
      close(78)

!
!    COEFFICIENTS FOR DIFFERENTIATION FOR NON-UNIFORM GRID
!
!    Q3 DIFFERENTIATION (CENTERED VARIABLE)
!

      am3ck(1)=0.d0
      ap3ck(1)=0.d0
      ac3ck(1)=1.d0
      am3ck(n3)=0.d0
      ap3ck(n3)=0.d0
      ac3ck(n3)=1.d0

      do kc=2,n3m
       km=kc-1
       kp=kc+1
       a33=dx3q/g3rc(kc)
       a33p=1.d0/g3rm(kc)
       a33m=1.d0/g3rm(km)
       ap3ck(kc)=a33*a33p
       am3ck(kc)=a33*a33m
       ac3ck(kc)=-(ap3ck(kc)+am3ck(kc))
      enddo

!
!    Q1/Q2 DIFFERENTIATION (STAGGERED VARIABLE)
!
!

      do kc=2,n3m-1
      kp=kc+1
      km=kc-1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      ap3sk(kc)=a33p
      am3sk(kc)=a33m
      ac3sk(kc)=-(ap3sk(kc)+am3sk(kc))
      enddo
!    
!    LOWER WALL BOUNDARY CONDITIONS (INSLWS SETS NO-SLIP vs STRESS-FREE WALL)
!    
      kc=1
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      ap3sk(kc)=a33p
      am3sk(kc)=0.d0
      ac3sk(kc)=-(a33p+inslws*a33m*2.d0)

!    
!    UPPER WALL BOUNDARY CONDITIONS (INSLWN SETS NO-SLIP vs STRESS-FREE WALL)
!    

      kc=n3m
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      am3sk(kc)=a33m
      ap3sk(kc)=0.d0
      ac3sk(kc)=-(a33m+inslwn*a33p*2.d0)

      am3ssk(1)=0.d0
      ap3ssk(1)=0.d0
      ac3ssk(1)=1.d0

!
!    TEMPERATURE DIFFERENTIATION (CENTERED VARIABLE)
!


      do kc=2,n3m
       kp=kc+1
       km=kc-1
       a33=dx3q/g3rc(kc)
       a33p=1.d0/g3rm(kc)
       a33m=1.d0/g3rm(km)
       ap3ssk(kc)=a33*a33p
       am3ssk(kc)=a33*a33m
       ac3ssk(kc)=-(ap3ssk(kc)+am3ssk(kc))
      enddo

      am3ssk(n3)=0.d0
      ap3ssk(n3)=0.d0
      ac3ssk(n3)=1.d0


      endif
      return                                                            
      end                                                               

