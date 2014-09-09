!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: StatReadReduceWrite.F90                        !
!    CONTAINS: subroutine StatReadReduceWrite             !
!                                                         ! 
!    PURPOSE: Reduce statistics across processors, read   !
!     old statistics if required and write them out       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine StatReadReduceWrite(var,filename,dsetname)
      use param
      use decomp_2d, only: xstart,xend,nrank
      use mpih
      implicit none
      character*30,intent(in) :: filename, dsetname
      real :: var(1:n3m,xstart(2):xend(2))
      real :: var_old(1:n3m)
      real :: var_new(1:n3m)



      call MPI_REDUCE(var,var_new,n3m, &
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

      if (nrank.eq.0) then

       if(starea.eq.1) then
        call HdfSerialReadReal1D(dsetname,filename,var_old,n3m)
        var_new = var_new + var_old
       endif

       call HdfSerialWriteReal1D(dsetname,filename,var_new,n3m)

      end if

       return
       end


