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
      use mpih
      implicit none
      character*30,intent(in) :: filename, dsetname
      real :: var(1:nxm)
      real :: var_old(1:nxm)



      call MpiSumReal1D(var,nxm)

      if (ismaster) then
       if(readstats) then
        call HdfSerialReadReal1D(dsetname,filename,var_old,nxm)
        var = var + var_old
       endif
       call HdfSerialWriteReal1D(dsetname,filename,var,nxm)
      end if

       return
       end


