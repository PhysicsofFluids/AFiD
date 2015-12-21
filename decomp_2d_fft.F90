!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This is the FFTW (version 3.x) implementation of the FFT library

module decomp_2d_fft
  
  use decomp_2d  ! 2D decomposition module
  
  implicit none

  INTEGER FFTW_ESTIMATE
  PARAMETER (FFTW_ESTIMATE=64)
  INTEGER FFTW_MEASURE
  PARAMETER (FFTW_MEASURE=0)
  INTEGER FFTW_FORWARD
  PARAMETER (FFTW_FORWARD=-1)
  INTEGER FFTW_BACKWARD
  PARAMETER (FFTW_BACKWARD=+1)
  
  private        ! Make everything private unless declared public

  ! engine-specific global variables
! integer, save :: plan_type = FFTW_EXHAUSTIVE
  integer, save :: plan_type = FFTW_ESTIMATE
!  integer, save :: plan_type = FFTW_MEASURE

  ! FFTW plans
  ! j=1,2,3 corresponds to the 1D FFTs in X,Y,Z direction, respectively
  ! For c2c transforms: 
  !     use plan(-1,j) for  forward transform; 
  !     use plan( 1,j) for backward transform;
  ! For r2c/c2r transforms:
  !     use plan(0,j) for r2c transforms;
  !     use plan(2,j) for c2r transforms;
  integer*8, save, public :: plan(-1:2,3)

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines
#include "fft_common.F90"

  ! Return a FFTW3 plan for multiple 1D c2c FFTs in X direction
  subroutine c2c_1m_x_plan(plan1, decomp, isign)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    complex(mytype), allocatable, dimension(:,:,:) :: a1

    allocate(a1(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))

    call dfftw_plan_many_dft(plan1, 1, decomp%xsz(1), &
         decomp%xsz(2)*decomp%xsz(3), a1, decomp%xsz(1), 1, &
         decomp%xsz(1), a1, decomp%xsz(1), 1, decomp%xsz(1), &
         isign, plan_type)

    deallocate(a1)

    return
  end subroutine c2c_1m_x_plan

  ! Return a FFTW3 plan for multiple 1D c2c FFTs in Y direction
  subroutine c2c_1m_y_plan(plan1, decomp, isign)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    complex(mytype), allocatable, dimension(:,:) :: a1

    ! Due to memory pattern of 3D arrays, 1D FFTs along Y have to be
    ! done one Z-plane at a time. So plan for 2D data sets here.

    allocate(a1(decomp%ysz(1),decomp%ysz(2)))

    call dfftw_plan_many_dft(plan1, 1, decomp%ysz(2), decomp%ysz(1), &
         a1, decomp%ysz(2), decomp%ysz(1), 1, a1, decomp%ysz(2), &
         decomp%ysz(1), 1, isign, plan_type)

    deallocate(a1)

    return
  end subroutine c2c_1m_y_plan


  ! Return a FFTW3 plan for multiple 1D c2c FFTs in Z direction
  subroutine c2c_1m_z_plan(plan1, decomp, isign)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    complex(mytype), allocatable, dimension(:,:,:) :: a1

    allocate(a1(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))

    call dfftw_plan_many_dft(plan1, 1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), a1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), 1, a1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), 1, isign, plan_type)

    deallocate(a1)

    return
  end subroutine c2c_1m_z_plan



  ! Return a FFTW3 plan for multiple 1D c2r FFTs in X direction
  subroutine c2r_1m_x_plan(plan1, decomp_sp, decomp_ph)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

    complex(mytype), allocatable, dimension(:,:,:) :: a1
    real(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_sp%xsz(1),decomp_sp%xsz(2),decomp_sp%xsz(3)))
    allocate(a2(decomp_ph%xsz(1),decomp_ph%xsz(2),decomp_ph%xsz(3)))
    call dfftw_plan_many_dft_c2r(plan1, 1, decomp_ph%xsz(1), &
         decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_sp%xsz(1), 1, &
         decomp_sp%xsz(1), a2, decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
         plan_type)
    deallocate(a1,a2)

    return
  end subroutine c2r_1m_x_plan


  ! Return a FFTW3 plan for multiple 1D r2c FFTs in Z direction
  subroutine r2c_1m_z_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    real(mytype), allocatable, dimension(:,:,:) :: a1
    complex(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_ph%zsz(1),decomp_ph%zsz(2),decomp_ph%zsz(3)))
    allocate(a2(decomp_sp%zsz(1),decomp_sp%zsz(2),decomp_sp%zsz(3)))
    call dfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, a2, decomp_sp%zsz(3), &
         decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, plan_type)
    deallocate(a1,a2)

    return
  end subroutine r2c_1m_z_plan


  ! Return a FFTW3 plan for multiple 1D c2r FFTs in Z direction
  subroutine c2r_1m_z_plan(plan1, decomp_sp, decomp_ph)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

    complex(mytype), allocatable, dimension(:,:,:) :: a1
    real(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_sp%zsz(1),decomp_sp%zsz(2),decomp_sp%zsz(3)))
    allocate(a2(decomp_ph%zsz(1),decomp_ph%zsz(2),decomp_ph%zsz(3)))

    call dfftw_plan_many_dft_c2r(plan1, 1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_sp%zsz(3), &
         decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, a2, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, plan_type)
    deallocate(a1,a2)     

    return
  end subroutine c2r_1m_z_plan

  ! Return a FFTW3 plan for multiple 1D r2c FFTs in X direction
  subroutine r2c_1m_x_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    real(mytype), allocatable, dimension(:,:,:) :: a1
    complex(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_ph%xsz(1),decomp_ph%xsz(2),decomp_ph%xsz(3)))
    allocate(a2(decomp_sp%xsz(1),decomp_sp%xsz(2),decomp_sp%xsz(3)))
    call dfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%xsz(1), &
         decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_ph%xsz(1), 1, &
         decomp_ph%xsz(1), a2, decomp_sp%xsz(1), 1, decomp_sp%xsz(1), &
         plan_type)
    deallocate(a1,a2)    

    return
  end subroutine r2c_1m_x_plan

  ! Return a FFTW3 plan for multiple 1D r2c FFTs in Y direction
  subroutine r2c_1m_y_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    real(mytype), allocatable, dimension(:,:,:) :: a1
    complex(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_ph%ysz(1),decomp_ph%ysz(2),decomp_ph%ysz(3)))
    allocate(a2(decomp_sp%ysz(1),decomp_sp%ysz(2),decomp_sp%ysz(3)))
    call dfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%ysz(1), &
         decomp_ph%ysz(2)*decomp_ph%ysz(3), a1, decomp_ph%ysz(1), 1, &
         decomp_ph%ysz(1), a2, decomp_sp%ysz(1), 1, decomp_sp%ysz(1), &
         plan_type)
    deallocate(a1,a2)    

    return
  end subroutine r2c_1m_y_plan

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine
!$  use omp_lib
    implicit none

!$  integer ierr

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the FFTW (version 3.x) engine *****'
       write(*,*) ' '
    end if

!$  call dfftw_init_threads(ierr)
!$ 
!$  call dfftw_plan_with_nthreads(omp_get_max_threads())

    if (format == PHYSICAL_IN_X) then

       ! For C2C transforms
!      call c2c_1m_x_plan(plan(-1,1), ph, FFTW_FORWARD )
!      call c2c_1m_y_plan(plan(-1,2), ph, FFTW_FORWARD )
!      call c2c_1m_z_plan(plan(-1,3), ph, FFTW_FORWARD )
!      call c2c_1m_z_plan(plan( 1,3), ph, FFTW_BACKWARD)
!      call c2c_1m_y_plan(plan( 1,2), ph, FFTW_BACKWARD)
!      call c2c_1m_x_plan(plan( 1,1), ph, FFTW_BACKWARD)
!      
!      ! For R2C/C2R tranforms
!      call r2c_1m_x_plan(plan(0,1), ph, sp)
!      call c2c_1m_y_plan(plan(0,2), sp, FFTW_FORWARD )
!      call c2c_1m_z_plan(plan(0,3), sp, FFTW_FORWARD )
!      call c2c_1m_z_plan(plan(2,3), sp, FFTW_BACKWARD)
!      call c2c_1m_y_plan(plan(2,2), sp, FFTW_BACKWARD)
!      call c2r_1m_x_plan(plan(2,1), sp, ph)

       
       ! For R2C/C2R tranforms
       call r2c_1m_y_plan(plan(1,1), ph, sp)
       call c2c_1m_y_plan(plan(1,2), sp, FFTW_BACKWARD)
       call c2r_1m_z_plan(plan(2,1), ph, sp)
       call c2c_1m_z_plan(plan(2,2), ph, FFTW_FORWARD)

    else if (format == PHYSICAL_IN_Z) then

       ! For C2C transforms
       call c2c_1m_z_plan(plan(-1,3), ph, FFTW_FORWARD )
       call c2c_1m_y_plan(plan(-1,2), ph, FFTW_FORWARD ) 
       call c2c_1m_x_plan(plan(-1,1), ph, FFTW_FORWARD )
       call c2c_1m_x_plan(plan( 1,1), ph, FFTW_BACKWARD)
       call c2c_1m_y_plan(plan( 1,2), ph, FFTW_BACKWARD)
       call c2c_1m_z_plan(plan( 1,3), ph, FFTW_BACKWARD)
       
       ! For R2C/C2R tranforms
       call r2c_1m_z_plan(plan(0,3), ph, sp)
       call c2c_1m_y_plan(plan(0,2), sp, FFTW_FORWARD )
       call c2c_1m_x_plan(plan(0,1), sp, FFTW_FORWARD )
       call c2c_1m_x_plan(plan(2,1), sp, FFTW_BACKWARD)
       call c2c_1m_y_plan(plan(2,2), sp, FFTW_BACKWARD)
       call c2r_1m_z_plan(plan(2,3), sp, ph)
       
    end if

    return
  end subroutine init_fft_engine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none

    integer :: i,j
    
    do j=1,3
       do i=-1,2
          call dfftw_destroy_plan(plan(i,j))
       end do
    end do

    return
  end subroutine finalize_fft_engine


  ! Following routines calculate multiple one-dimensional FFTs to form 
  ! the basis of three-dimensional FFTs.

  ! c2c transform, multiple 1D FFTs in x direction
  subroutine c2c_1m_x(inout, plan1)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer*8, intent(IN) :: plan1

    call dfftw_execute_dft(plan1, inout, inout)

    return
  end subroutine c2c_1m_x


  ! c2c transform, multiple 1D FFTs in y direction
  subroutine c2c_1m_y(inout, plan1)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer*8, intent(IN) :: plan1

    integer :: k, s3

    ! transform on one Z-plane at a time
    s3 = size(inout,3)
    do k=1,s3
       call dfftw_execute_dft(plan1, inout(:,:,k), inout(:,:,k))
    end do

    return
  end subroutine c2c_1m_y

  ! c2c transform, multiple 1D FFTs in z direction
  subroutine c2c_1m_z(inout, plan1)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer*8, intent(IN) :: plan1

       call dfftw_execute_dft(plan1, inout, inout)

    return
  end subroutine c2c_1m_z

  ! r2c transform, multiple 1D FFTs in x direction
  subroutine r2c_1m_x(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    call dfftw_execute_dft_r2c(plan(0,1), input, output)

    return

  end subroutine r2c_1m_x

  ! r2c transform, multiple 1D FFTs in x direction
  subroutine r2c_1m_y(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    call dfftw_execute_dft_r2c(plan(1,1), input, output)

    return

  end subroutine r2c_1m_y

  ! r2c transform, multiple 1D FFTs in z direction
  subroutine r2c_1m_z(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    call dfftw_execute_dft_r2c(plan(0,3), input, output)

    return

  end subroutine r2c_1m_z

  ! c2r transform, multiple 1D FFTs in x direction
  subroutine c2r_1m_x(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    call dfftw_execute_dft_c2r(plan(2,1), input, output)

    return

  end subroutine c2r_1m_x

  ! c2r transform, multiple 1D FFTs in z direction
  subroutine c2r_1m_z(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN) :: input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    call dfftw_execute_dft_c2r(plan(2,1), input, output)

    return

  end subroutine c2r_1m_z
end module decomp_2d_fft
