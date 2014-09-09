      subroutine quit(tin,cond)
      use hdf5
      use mpih
      use param
      use decomp_2d, only: nrank
      use decomp_2d_fft
      implicit none
      integer, intent(in) :: cond
      integer :: hdf_error
      real :: tin(3)

      tin(3) = MPI_WTIME()
      if(nrank.eq.0) then
          write(6,*) 'Total Iteration Time = ',tin(3) -tin(2),' sec.'
      endif

      if(cond.eq.1) then
        if (statcal) call ststwr
        call continua
      endif

      call closefi
      call DeallocateVariables
      call h5close_f(hdf_error)
      call decomp_2d_fft_finalize

      stop

      end subroutine quit

