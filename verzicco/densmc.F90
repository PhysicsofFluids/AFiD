
!***********************************************************************

!***********************************************************************
      subroutine densmc
      use param
      use local_arrays, only: dens
      use mpih
      use decomp_2d, only: nrank,xstart,xend,xsize,DECOMP_2D_COMM_CART_X
      implicit none
      integer :: j,i
      real :: grtlow,grtupp
      real ::  my_anussupp,my_anusslow
      real ::  anussupp,anusslow
      real :: del,deln,my_surface
  
!
!     COMPUTATION OF THE NUSSELT NUMBER AT THE 
!     LOWER AND UPPER WALLS
!

      anusslow = 0.d0
      my_anusslow = 0.d0
      del  = 1.0/(zz(2)-zz(1))
      anussupp = 0.d0
      my_anussupp = 0.d0
      deln = 1.0/(zz(n3)-zz(n3m))

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,dens,del,deln) &
!$OMP   SHARED(n3m,n3) &
!$OMP   PRIVATE(i,j,grtlow,grtupp) &
!$OMP   REDUCTION(+:my_anusslow) &
!$OMP   REDUCTION(+:my_anussupp)
      do i=xstart(3),xend(3)
         do j=xstart(2),xend(2)
           grtlow = (dens(1,j,i)-dens(2,j,i))*del
           my_anusslow = my_anusslow + grtlow
           grtupp = (dens(n3m,j,i)-dens(n3,j,i))*deln
           my_anussupp = my_anussupp + grtupp
        enddo
      end do
!$OMP END PARALLEL DO

      my_surface = 1.0/(n1m*n2m)
      my_anusslow = my_anusslow * my_surface
      my_anussupp = my_anussupp * my_surface

      call MPI_REDUCE(my_anusslow,anusslow,1,MDP,MPI_SUM,0,DECOMP_2D_COMM_CART_X,ierr)
      call MPI_REDUCE(my_anussupp,anussupp,1,MDP,MPI_SUM,0,DECOMP_2D_COMM_CART_X,ierr)

      if(nrank.eq.0) then
       write(97,546) time, anusslow, anussupp
 546   format(4(1x,e14.6))
      endif

      return         
      end                                                               
