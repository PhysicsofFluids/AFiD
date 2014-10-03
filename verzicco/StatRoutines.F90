!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: stst.F90                                       !
!    CONTAINS: subroutine CalcStats,WriteStats            !
!                                                         ! 
!    PURPOSE: Calculates and writes out statistics for    !
!     the flow field. All quantities are averaged in the  !
!     two horizontal (homogeneous) directions             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcStats
      use param
      use local_arrays, only: q1,q2,q3,temp
      use decomp_2d, only: xstart,xend
      use stat_arrays
      use mpih
      implicit none
      real :: usn1m,usn2m
      integer :: i,j,k

      timeint_cdsp = timeint_cdsp + 1

      usn1m = 1.0/nzm
      usn2m = 1.0/nym

      do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
          do k=1,nxm
               q1_me(k) = q1_me(k) + q1(k,j,i)*usn1m*usn2m
               q2_me(k) = q2_me(k) + q2(k,j,i)*usn1m*usn2m
               q3_me(k) = q3_me(k) + q3(k,j,i)*usn1m*usn2m
               temp_me(k) = temp_me(k) + temp(k,j,i)*usn1m*usn2m
               q1_rms(k) = q1_rms(k) + q1(k,j,i)**2*usn1m*usn2m
               q2_rms(k) = q2_rms(k) + q2(k,j,i)**2*usn1m*usn2m
               q3_rms(k) = q3_rms(k) + q3(k,j,i)**2*usn1m*usn2m
               temp_rms(k) = temp_rms(k) +  &
     &                               temp(k,j,i)**2*usn1m*usn2m
               tempq3_me(k) = tempq3_me(k) +  &
     &                         temp(k,j,i)*q3(k,j,i)*usn1m*usn2m
            end do
         end do
      end do

      return  
      end
!    
!***********************************************************************
      subroutine WriteStats
      use mpih
      use param
      use stat_arrays
      use hdf5

      implicit none

      integer :: timeint_cdsp_old

      character*30 dsetname_q1me
      character*30 dsetname_q2me
      character*30 dsetname_q3me
      character*30 dsetname_tempme
      character*30 dsetname_q1rms
      character*30 dsetname_q2rms
      character*30 dsetname_q3rms
      character*30 dsetname_temprms
      character*30 dsetname_tempq3me
      character*30 dsetname_dissth
      character*30 dsetname_disste
      character*30 filnam,dsetname
      logical :: fexist

      filnam = trim('stafield_master.h5')

      dsetname_q1me = trim('q1_mean')
      dsetname_q2me = trim('q2_mean')
      dsetname_q3me = trim('q3_mean')
      dsetname_tempme = trim('temp_mean.h5')
      dsetname_q1rms = trim('q1_rms.h5')
      dsetname_q2rms = trim('q2_rms.h5')
      dsetname_q3rms = trim('q3_rms.h5')
      dsetname_temprms = trim('temp_rms.h5')
      dsetname_tempq3me = trim('tempq3_mean.h5')
      dsetname_dissth = trim('dissth.h5')
      dsetname_disste = trim('disste.h5')

      dsetname = trim('averaging_time')

      inquire(file=filnam,exist=fexist)
      if (.not.fexist) then 
        if(ismaster) write(6,*) 'Unable to read statistical files'
        if(ismaster) write(6,*) 'Restarting statistics from zero' 
        readstats=.false.
      end if
       

      if (ismaster) then
       if(readstats) then
        call HdfSerialReadIntScalar(dsetname,filnam,timeint_cdsp_old)
        timeint_cdsp = timeint_cdsp + timeint_cdsp_old
       else 
        call HdfCreateBlankFile(filnam)
       endif
      end if

      call StatReadReduceWrite(q1_me,filnam,dsetname_q1me)
      call StatReadReduceWrite(q2_me,filnam,dsetname_q2me)
      call StatReadReduceWrite(q3_me,filnam,dsetname_q3me)
      call StatReadReduceWrite(temp_me,filnam,dsetname_tempme)

      call StatReadReduceWrite(q1_rms,filnam,dsetname_q1rms)
      call StatReadReduceWrite(q2_rms,filnam,dsetname_q2rms)
      call StatReadReduceWrite(q3_rms,filnam,dsetname_q3rms)
      call StatReadReduceWrite(temp_rms,filnam,dsetname_temprms)
 
      call StatReadReduceWrite(tempq3_me,filnam,dsetname_tempq3me)

      if(disscal) then 
       call StatReadReduceWrite(dissth,filnam,dsetname_dissth)
       call StatReadReduceWrite(disste,filnam,dsetname_disste)
      end if

      if (ismaster) then

       call HdfSerialWriteIntScalar(dsetname,filnam,timeint_cdsp)

       dsetname = trim('Rayleigh Number')
       call HdfSerialWriteRealScalar(dsetname,filnam,ray)

       dsetname = trim('Prandtl Number')
       call HdfSerialWriteRealScalar(dsetname,filnam,pra)

       dsetname = trim('X_cordin')
       call HdfSerialWriteReal1D(dsetname,filnam,zm,1,nxm)

      endif

      return  
      end
! 
