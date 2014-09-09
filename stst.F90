!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: stst.F90                                       !
!    CONTAINS: subroutine stst,ststwr                     !
!                                                         ! 
!    PURPOSE: Calculates and writes out statistics for    !
!     the flow field. All quantities are averaged in the  !
!     two horizontal (homogeneous) directions             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine stst
      use param
      use local_arrays, only: q1,q2,q3,dens
      use decomp_2d, only: xstart,xend
      use stat_arrays
      use mpih
      implicit none
      real :: usn1m,usn2m
      integer :: i,j,k

      timeint_cdsp = timeint_cdsp + 1

      usn1m = 1.0/n1m
      usn2m = 1.0/n2m

      do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
          do k=1,n3m
               q1_me(k) = q1_me(k) + q1(k,j,i)*usn1m*usn2m
               q2_me(k) = q2_me(k) + q2(k,j,i)*usn1m*usn2m
               q3_me(k) = q3_me(k) + q3(k,j,i)*usn1m*usn2m
               dens_me(k) = dens_me(k) + dens(k,j,i)*usn1m*usn2m
               q1_rms(k) = q1_rms(k) + q1(k,j,i)**2*usn1m*usn2m
               q2_rms(k) = q2_rms(k) + q2(k,j,i)**2*usn1m*usn2m
               q3_rms(k) = q3_rms(k) + q3(k,j,i)**2*usn1m*usn2m
               dens_rms(k) = dens_rms(k) +  &
     &                               dens(k,j,i)**2*usn1m*usn2m
               densq3_me(k) = densq3_me(k) +  &
     &                         dens(k,j,i)*q3(k,j,i)*usn1m*usn2m
            end do
         end do
      end do

      return  
      end
!    
!***********************************************************************
      subroutine ststwr
      use mpih
      use param
      use stat_arrays
      use hdf5
      use decomp_2d, only: nrank

      implicit none

      integer hdf_error

      integer(HID_T) :: file_id

      integer(HSIZE_T) :: dims_grid(1)
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HID_T) :: plist_full

      integer :: ndims,timeint_cdsp_old

      character*30 filename_q1me,dsetname_q1me
      character*30 filename_q2me,dsetname_q2me
      character*30 filename_q3me,dsetname_q3me
      character*30 filename_densme,dsetname_densme
      character*30 filename_q1rms,dsetname_q1rms
      character*30 filename_q2rms,dsetname_q2rms
      character*30 filename_q3rms,dsetname_q3rms
      character*30 filename_densrms,dsetname_densrms
      character*30 filename_densq3me,dsetname_densq3me
      character*30 filename_dissth,dsetname_dissth
      character*30 filename_disste,dsetname_disste
      character*30 filnamgrid

      filename_q1me = trim('stats/q1_mean.h5')
      filename_q2me = trim('stats/q2_mean.h5')
      filename_q3me = trim('stats/q3_mean.h5')
      filename_densme = trim('stats/dens_mean.h5')
      filename_q1rms = trim('stats/q1_rms.h5')
      filename_q2rms = trim('stats/q2_rms.h5')
      filename_q3rms = trim('stats/q3_rms.h5')
      filename_densrms = trim('stats/dens_rms.h5')
      filename_densq3me = trim('stats/densq3_mean.h5')
      filename_dissth = trim('stats/dissth.h5')
      filename_disste = trim('stats/disste.h5')

      dsetname_q1me = trim('q1_mean')
      dsetname_q2me = trim('q2_mean')
      dsetname_q3me = trim('q3_mean')
      dsetname_densme = trim('dens_mean.h5')
      dsetname_q1rms = trim('q1_rms.h5')
      dsetname_q2rms = trim('q2_rms.h5')
      dsetname_q3rms = trim('q3_rms.h5')
      dsetname_densrms = trim('dens_rms.h5')
      dsetname_densq3me = trim('densq3_mean.h5')
      dsetname_dissth = trim('dissth.h5')
      dsetname_disste = trim('disste.h5')

      call StatReadReduceWrite(q1_me,filename_q1me,dsetname_q1me)
      call StatReadReduceWrite(q2_me,filename_q2me,dsetname_q2me)
      call StatReadReduceWrite(q3_me,filename_q3me,dsetname_q3me)
      call StatReadReduceWrite(dens_me,filename_densme,dsetname_densme)

      call StatReadReduceWrite(q1_rms,filename_q1rms,dsetname_q1rms)
      call StatReadReduceWrite(q2_rms,filename_q2rms,dsetname_q2rms)
      call StatReadReduceWrite(q3_rms,filename_q3rms,dsetname_q3rms)
      call StatReadReduceWrite(dens_rms,filename_densrms,dsetname_densrms)
 
      call StatReadReduceWrite(densq3_me,filename_densq3me,dsetname_densq3me)

      if(balcal) then 
       call StatReadReduceWrite(dissth,filename_dissth,dsetname_dissth)
       call StatReadReduceWrite(disste,filename_disste,dsetname_disste)
      end if


      filnamgrid = 'stafield_master.h5'
      ndims=1

      if (nrank.eq.0) then
       if(starea.eq.1) then
        dims_grid(1)=1
  
        call h5fopen_f(filnamgrid, H5F_ACC_RDONLY_F, file_id, hdf_error)
  
        call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
  
        call h5dopen_f(file_id, 'averaging_time', dset_grid, hdf_error)
  
        call h5dread_f(dset_grid, H5T_NATIVE_INTEGER, timeint_cdsp_old, &
               dims_grid,hdf_error)
  
        call h5dclose_f(dset_grid, hdf_error)
        call h5sclose_f(dspace_grid, hdf_error)
  
        call h5pclose_f(plist_full, hdf_error)
        call h5fclose_f(file_id, hdf_error)
  
        timeint_cdsp = timeint_cdsp + timeint_cdsp_old

       endif

!RO   Write the grid & statistics information
!RO   only if master process


      call h5fcreate_f(filnamgrid,H5F_ACC_TRUNC_F, file_id, hdf_error)

!RO   Write amount of averages 

      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'averaging_time', H5T_NATIVE_INTEGER, &
     &                dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_INTEGER, timeint_cdsp, &
     &       dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!RO   Write Reynolds number

      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'Ra', H5T_NATIVE_DOUBLE, &
     &                dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ray, &
     &       dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!EP   Write Prandtl number

      dims_grid(1)=1
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_id, 'Pr', H5T_NATIVE_DOUBLE, &
     &                dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, pra, &
     &       dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!RO   Write the grid information 

      dims_grid(1)=n3m
      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_id, 'Z_cordin', H5T_NATIVE_DOUBLE, &
     &                dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m), &
     &        dims_grid, hdf_error)


      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

!RO   Close file

      call h5fclose_f(file_id, hdf_error)

      endif

      return  
      end
! 
