!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: StatRoutines.F90                               !
!    CONTAINS: subroutine CalcStats,WriteStats            !
!                                                         ! 
!    PURPOSE: Calculates and writes out statistics for    !
!     the flow field. All quantities are averaged in the  !
!     two horizontal (homogeneous) directions             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcStats
      use param
      use local_arrays, only: vz,vy,vx,temp
      use decomp_2d, only: xstart,xend
      use stat_arrays
      use mpih
      implicit none
      real :: usnzm,usnym
      integer :: i,j,k

      nstatsamples = nstatsamples + 1

      usnym = 1.0/nym
      usnzm = 1.0/nzm

      do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
          do k=1,nxm
               vx_me(k) = vx_me(k) + vx(k,j,i)*usnzm*usnym
               vy_me(k) = vy_me(k) + vy(k,j,i)*usnzm*usnym
               vz_me(k) = vz_me(k) + vz(k,j,i)*usnzm*usnym
               temp_me(k) = temp_me(k) + temp(k,j,i)*usnzm*usnym
               vx_rms(k) = vx_rms(k) + vx(k,j,i)**2*usnzm*usnym
               vy_rms(k) = vy_rms(k) + vy(k,j,i)**2*usnzm*usnym
               vz_rms(k) = vz_rms(k) + vz(k,j,i)**2*usnzm*usnym
               temp_rms(k) = temp_rms(k) +  &
     &                               temp(k,j,i)**2*usnzm*usnym
               tempvx_me(k) = tempvx_me(k) +  &
     &                         temp(k,j,i)*vx(k,j,i)*usnzm*usnym
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

      integer :: nstatsamples_old

      character*30 dsetname_vxme
      character*30 dsetname_vyme
      character*30 dsetname_vzme
      character*30 dsetname_tempme

      character*30 dsetname_vxrms
      character*30 dsetname_vyrms
      character*30 dsetname_vzrms
      character*30 dsetname_temprms

      character*30 dsetname_tempvxme
      character*30 dsetname_dissth
      character*30 dsetname_disste
      character*30 filnam,dsetname
      logical :: fexist

      filnam = trim('stafield_master.h5')

      dsetname_vxme = trim('vx_mean')
      dsetname_vyme = trim('vy_mean')
      dsetname_vzme = trim('vz_mean')
      dsetname_tempme = trim('temp_mean.h5')

      dsetname_vxrms = trim('vx_rms.h5')
      dsetname_vyrms = trim('vy_rms.h5')
      dsetname_vzrms = trim('vz_rms.h5')
      dsetname_temprms = trim('temp_rms.h5')

      dsetname_tempvxme = trim('tempvx_mean.h5')
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
        call HdfSerialReadIntScalar(dsetname,filnam,nstatsamples_old)
        nstatsamples = nstatsamples + nstatsamples_old
       else 
        call HdfCreateBlankFile(filnam)
       endif
      end if

      call StatReadReduceWrite(vx_me,filnam,dsetname_vxme)
      call StatReadReduceWrite(vy_me,filnam,dsetname_vyme)
      call StatReadReduceWrite(vz_me,filnam,dsetname_vzme)
      call StatReadReduceWrite(temp_me,filnam,dsetname_tempme)

      call StatReadReduceWrite(vx_rms,filnam,dsetname_vxrms)
      call StatReadReduceWrite(vy_rms,filnam,dsetname_vyrms)
      call StatReadReduceWrite(vz_rms,filnam,dsetname_vzrms)
      call StatReadReduceWrite(temp_rms,filnam,dsetname_temprms)
 
      call StatReadReduceWrite(tempvx_me,filnam,dsetname_tempvxme)

      if(disscal) then 
       call StatReadReduceWrite(dissth,filnam,dsetname_dissth)
       call StatReadReduceWrite(disste,filnam,dsetname_disste)
      end if

      if (ismaster) then

       call HdfSerialWriteIntScalar(dsetname,filnam,nstatsamples)

       dsetname = trim('X_cordin')
       call HdfSerialWriteReal1D(dsetname,filnam,xm,1,nxm)

       dsetname = trim('Rayleigh Number')
       call HdfSerialWriteRealScalar(dsetname,filnam,ray)

       dsetname = trim('Prandtl Number')
       call HdfSerialWriteRealScalar(dsetname,filnam,pra)


      endif

      return  
      end
! 
