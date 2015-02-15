!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: HdfRoutines.F90                                !
!    CONTAINS: subroutines MPI*                           !
!                                                         ! 
!    PURPOSE: Wrappers for MPI Routines                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine MpiBcastInt(n)
      use mpih
      implicit none
      integer, intent(in) :: n
      
      call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      return
      end subroutine MpiBcastInt

!==============================================================================

      subroutine MpiBcastReal(n)
      use mpih
      implicit none
      real, intent(in) :: n
      
      call MPI_BCAST(n,1,MDP,0,MPI_COMM_WORLD,ierr)

      return
      end subroutine MpiBcastReal
!==============================================================================

      subroutine MpiBarrier
      use mpih
      implicit none
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end subroutine MpiBarrier

!==============================================================================

      subroutine MpiSumRealScalar(var)
      use mpih
      implicit none
      real, intent(inout) :: var
      real :: buf
      
       call MPI_REDUCE(var,buf,1, &
        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

       var = buf

      return
      end subroutine MpiSumRealScalar
!==============================================================================

      subroutine MpiMaxRealScalar(var)
      use mpih
      implicit none
      real, intent(inout) :: var
      real :: buf
      
       call MPI_REDUCE(var,buf,1, &
        MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
 
       var = buf

      return
      end subroutine MpiMaxRealScalar
!==============================================================================

      subroutine MpiMinRealScalar(var)
      use mpih
      implicit none
      real, intent(inout) :: var
      real :: buf
      
       call MPI_REDUCE(var,buf,1, &
        MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
 
       var = buf

      return
      end subroutine MpiMinRealScalar
!==============================================================================

      subroutine MpiSumReal1D(var,sz)
      use mpih
      implicit none
      integer, intent(in) :: sz
      real, intent(inout), dimension(1:sz) :: var
      real, dimension(1:sz) :: buf
      
       call MPI_REDUCE(var,buf,sz, &
        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 
       var = buf

      return
      end subroutine MpiSumReal1D
!==============================================================================

      subroutine MpiAbort
      use mpih
      implicit none
      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)

      return
      end subroutine MpiAbort
