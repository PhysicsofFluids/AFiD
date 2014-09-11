!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitPressureSolver.F90                         !
!    CONTAINS: subroutine InitPressureSolver              !
!                                                         ! 
!    PURPOSE: Initialization routines. Compute the metric !
!     terms and modified wavenumbers for the pressure     !
!     correction                                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine InitPressureSolver
      use param
      use decomp_2d_fft
      implicit none
      integer :: n2mh,n2mp,j,i,n1mh,n1mp
      integer  :: kc,km,kp
      real :: ugmmm,a33icc,a33icp
    
!    Initialize wave number definitions

      n1mh=n1m/2+1
      n1mp=n1mh+1
      n2mh=n2m/2+1
      n2mp=n2mh+1

      do i=1,n1mh
        ao(i)=(i-1)*2.d0*pi
      enddo
      do i=n1mp,n1m
        ao(i)=-(n1m-i+1)*2.d0*pi
      enddo
      do i=1,n1m
        ak1(i)=2.d0*(1.d0-dcos(ao(i)/n1m))*(float(n1m)/rext)**2
      enddo

      do j=1,n2mh
        ap(j)=(j-1)*2.d0*pi
      enddo
      do j=n2mp,n2m
        ap(j)=-(n2m-j+1)*2.d0*pi
      enddo
      do j=1,n2m
        ak2(j)=2.d0*(1.d0-dcos(ap(j)/n2m))*(float(n2m)/rext2)**2
      enddo

!RO   Initialize Tridiagonal matrices for Poisson solver
      do kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        a33icc=kmc(kc)*dx3q/g3rc(kc)
        a33icp=kpc(kc)*dx3q/g3rc(kp)
        ugmmm=1.0d0/g3rm(kc)
        amphk(kc)=a33icc*ugmmm
        apphk(kc)=a33icp*ugmmm
        acphk(kc)=-(amphk(kc)+apphk(kc))
      enddo

!EP   Initialize planning
      call decomp_2d_fft_init

      return
      end
      
