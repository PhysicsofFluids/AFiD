!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
! Increased efficiency and adapted to 3DRB by Erwin P. van der Poel 2013
!=======================================================================

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Halo cell support for neighbouring pencils to exchange data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_halo_real(in, level)

    implicit none

    integer, intent(IN) :: level      ! levels of halo cells required
    real(mytype), dimension(:,:,:), intent(INOUT) :: in    

    ! starting/ending index of array with halo cells
    integer :: xs, ys, zs, xe, ye, ze

    integer :: s1, s2, s3, ierror
    integer :: data_type

    integer :: icount, ilength, ijump 
    integer :: halo12             
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_n, tag_s, tag_t, tag_b

#ifdef HALO_DEBUG
      real :: cksum,cksum2,cksum3,cksum4
#endif

    data_type = real_type

#include "halo_common.f90"

    return
  end subroutine update_halo_real


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! To support halo-cell exchange:
  !   find the MPI ranks of neighbouring pencils
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_neighbour

    integer :: ierror

    ! For X-pencil
    neighbour(1,1) = MPI_PROC_NULL               ! east
    neighbour(1,2) = MPI_PROC_NULL               ! west
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, &
         neighbour(1,4), neighbour(1,3), ierror) ! north & south
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, &
         neighbour(1,6), neighbour(1,5), ierror) ! top & bottom

    ! For Y-pencil
!   call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 0, 1, &
!        neighbour(2,2), neighbour(2,1), ierror) ! east & west
!   neighbour(2,3) = MPI_PROC_NULL               ! north
!   neighbour(2,4) = MPI_PROC_NULL               ! south
!   call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 1, 1, &
!        neighbour(2,6), neighbour(2,5), ierror) ! top & bottom

!   ! For Z-pencil
!   call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 0, 1, &
!        neighbour(3,2), neighbour(3,1), ierror) ! east & west
!   call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 1, 1, &
!        neighbour(3,4), neighbour(3,3), ierror) ! north & south
!   neighbour(3,5) = MPI_PROC_NULL               ! top
!   neighbour(3,6) = MPI_PROC_NULL               ! bottom

    return
  end subroutine init_neighbour
