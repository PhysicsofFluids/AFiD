!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: QuitRoutine.F90                                !
!    CONTAINS: subroutine QuitRoutine, NotifyError        !
!                                                         ! 
!    PURPOSE: Routines to exit the program and write the  !
!     data if necessary                                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine QuitRoutine(tin,normalexit,errorcode)
      use hdf5
      use mpih
      use param
      use decomp_2d, only: nrank, decomp_2d_finalize
      use decomp_2d_fft
      implicit none
      logical, intent(in) :: normalexit
      integer :: errorcode
      real :: tin(3)

      if(errorcode.ne.100) then !EP skip if already finalized

      tin(3) = MPI_WTIME()
      if(ismaster) then
       call NotifyError(errorcode) 
      endif

      if(normalexit) then
        if(nrank.eq.0) write(6,'(a,f10.2,a)') 'Total Iteration Time = ',tin(3) -tin(2),' sec.'
        if (statcal) call WriteStats
        call WriteFlowField
      else
        call MPI_Abort(MPI_COMM_WORLD,1)
      endif

      call DeallocateVariables
      call HdfClose
      call decomp_2d_fft_finalize
      call decomp_2d_finalize

      endif

      end subroutine QuitRoutine


      subroutine NotifyError(errorcode)
      use param
      implicit none
      integer, intent(in) :: errorcode

      if(errorcode.eq.166) then 
        write(6,168) dt 
168     format(10x,'dt too small, DT= ',e14.7)
      else if(errorcode.eq.165) then
        write(6,164) 
164     format(10x,'cfl too large  ')
      else if(errorcode.eq.266) then
        write(6,268)
268     format(10x,'velocities diverged')
      else if(errorcode.eq.169) then
        write(6,178) 
178     format(10x,'too large local residue for mass conservation at:')
        call LocateLargeDivergence
      else if(errorcode.eq.333) then
         write(*,*) "time greater than tmax"
         write(*,*) "statistics and continuation updated"
      else if(errorcode.eq.334) then
         write(*,*) "walltime greater than walltimemax"
         write(*,*) "statistics and continuation updated"
      else 
         write(*,*) "Maximum number of timesteps reached"
      end if

      return
      end

     
