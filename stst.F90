!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: stst.F90                                       !
!    CONTAINS: subroutine stst,initstst,ststwr            !
!                                                         ! 
!    PURPOSE: Calculates and writes out statistics for    !
!     the flow field. All quantities are averaged in the  !
!     two horizontal (homogeneous) directions             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine stst
#ifdef STATS
      use param
      use local_arrays, only: q1,q2,q3,dens
      use decomp_2d, only: xstart,xend,nrank
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

#endif
      return  
      end
!    
!***********************************************************************
      subroutine ststwr
#ifdef STATS
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
      character*30 filnamgrid

      real,dimension(1:n3m) :: my_q1_me,my_q2_me,my_q3_me,my_dens_me
      real,dimension(1:n3m) :: my_q1_rms,my_q2_rms,my_q3_rms
      real,dimension(1:n3m) :: my_dens_rms,my_densq3_me

      real,dimension(1:n3m) :: q1_me_old,q2_me_old,q3_me_old,dens_me_old
      real,dimension(1:n3m) :: q1_rms_old,q2_rms_old,q3_rms_old
      real,dimension(1:n3m) :: dens_rms_old,densq3_me_old
#ifdef BALANCE
      real,dimension(1:n3m) :: dissth_old,disste_old
      real,dimension(1:n3m) :: my_dissth,my_disste
      character*30 filename_dissth,dsetname_dissth
      character*30 filename_disste,dsetname_disste
#endif
      

!EP   REDUCE

      my_q1_me=q1_me
      my_q2_me=q2_me
      my_q3_me=q3_me

      my_q1_rms=q1_rms
      my_q2_rms=q2_rms
      my_q3_rms=q3_rms

      my_dens_me=dens_me
      my_dens_rms=dens_rms

      my_densq3_me=densq3_me


      call MPI_REDUCE(my_q1_me,q1_me,n3m,MPI_DOUBLE_PRECISION,MPI_SUM, &
     & 0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q2_me,q2_me,n3m,MPI_DOUBLE_PRECISION,MPI_SUM, &
     & 0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q3_me,q3_me,n3m,MPI_DOUBLE_PRECISION,MPI_SUM, &
     & 0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_dens_me,dens_me,n3m, &
        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q1_rms,q1_rms,n3m, &
        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q2_rms,q2_rms,n3m, &
       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q3_rms,q3_rms, &
        n3m,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_dens_rms,dens_rms,n3m,MPI_DOUBLE_PRECISION, &
     & MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_densq3_me,densq3_me,n3m,MPI_DOUBLE_PRECISION, &
     & MPI_SUM,0,MPI_COMM_WORLD,ierr)
  
#ifdef BALANCE
      my_dissth=dissth
      my_disste=disste

      call MPI_REDUCE(my_dissth,dissth,n3m,MPI_DOUBLE_PRECISION, &
     & MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_disste,disste,n3m,MPI_DOUBLE_PRECISION, &
     & MPI_SUM,0,MPI_COMM_WORLD,ierr)
#endif

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

      if(starea.eq.1) then
      call hdf_read_serial_1d(dsetname_q1me,filename_q1me,n3m,q1_me_old)
      call hdf_read_serial_1d(dsetname_q2me,filename_q2me,n3m,q2_me_old)
      call hdf_read_serial_1d(dsetname_q3me,filename_q3me,n3m,q3_me_old)
      call hdf_read_serial_1d(dsetname_densme,filename_densme, &
                              n3m,dens_me_old)
      call hdf_read_serial_1d(dsetname_q1rms,filename_q1rms, &
     &                        n3m,q1_rms_old)
      call hdf_read_serial_1d(dsetname_q2rms,filename_q2rms, &
     &                        n3m,q2_rms_old)
      call hdf_read_serial_1d(dsetname_q3rms,filename_q3rms, &
     &                        n3m,q3_rms_old)
      call hdf_read_serial_1d(dsetname_densrms,filename_densrms, &
                              n3m,dens_rms_old)
      call hdf_read_serial_1d(dsetname_densq3me,filename_densq3me, &
                              n3m,densq3_me_old)
#ifdef BALANCE
      call hdf_read_serial_1d(dsetname_dissth,filename_dissth, &
     &                        n3m,dissth_old)
      call hdf_read_serial_1d(dsetname_disste,filename_disste, &
     &                        n3m,disste_old)
#endif

      ndims=1
      dims_grid(1)=1

      filnamgrid = 'stafield_master.h5'
      call h5fopen_f(filnamgrid, H5F_ACC_RDONLY_F, file_id, hdf_error)

      call h5screate_simple_f(ndims, dims_grid, dspace_grid, hdf_error)

      call h5dopen_f(file_id, 'averaging_time', dset_grid, hdf_error)

      call h5dread_f(dset_grid, H5T_NATIVE_INTEGER, timeint_cdsp_old, &
     &       dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5pclose_f(plist_full, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      q1_me = q1_me * timeint_cdsp + q1_me_old * timeint_cdsp_old
      q2_me = q2_me * timeint_cdsp + q2_me_old * timeint_cdsp_old
      q3_me = q3_me * timeint_cdsp + q3_me_old * timeint_cdsp_old
      dens_me = dens_me * timeint_cdsp + dens_me_old * timeint_cdsp_old
      q1_rms = q1_rms * timeint_cdsp + q1_rms_old * timeint_cdsp_old
      q2_rms = q2_rms * timeint_cdsp + q2_rms_old * timeint_cdsp_old
      q3_rms = q3_rms * timeint_cdsp + q3_rms_old * timeint_cdsp_old
      dens_rms = dens_rms * timeint_cdsp +  &
     &                         dens_rms_old * timeint_cdsp_old
      densq3_me = densq3_me * timeint_cdsp + &
     &                         densq3_me_old * timeint_cdsp_old
#ifdef BALANCE
      dissth = dissth * timeint_cdsp + dissth_old * timeint_cdsp_old
      disste = disste * timeint_cdsp + disste_old * timeint_cdsp_old
#endif
      
      timeint_cdsp = timeint_cdsp + timeint_cdsp_old

      endif

      if (nrank.eq.0) then

      call hdf_write_serial_1d(dsetname_q1me,filename_q1me,n3m,q1_me)
      call hdf_write_serial_1d(dsetname_q2me,filename_q2me,n3m,q2_me)
      call hdf_write_serial_1d(dsetname_q3me,filename_q3me,n3m,q3_me)
      call hdf_write_serial_1d(dsetname_densme,filename_densme, &
                              n3m,dens_me)
      call hdf_write_serial_1d(dsetname_q1rms,filename_q1rms, &
     &                        n3m,q1_rms)
      call hdf_write_serial_1d(dsetname_q2rms,filename_q2rms, &
     &                        n3m,q2_rms)
      call hdf_write_serial_1d(dsetname_q3rms,filename_q3rms, &
     &                        n3m,q3_rms)
      call hdf_write_serial_1d(dsetname_densrms,filename_densrms, &
                              n3m,dens_rms)
      call hdf_write_serial_1d(dsetname_densq3me,filename_densq3me, &
                              n3m,densq3_me)
#ifdef BALANCE
      call hdf_write_serial_1d(dsetname_dissth,filename_dissth, &
     &                        n3m,dissth)
      call hdf_write_serial_1d(dsetname_disste,filename_disste, &
     &                        n3m,disste)
#endif

!RO   Write the grid & statistics information
!RO   only if master process

      ndims=1

      filnamgrid = 'stafield_master.h5'
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

#endif
      return  
      end
! 
!***********************************************************************

      subroutine initstst
#ifdef STATS
      use param
      use stat_arrays
      implicit none
      integer :: k

!EP   Read or initialize stat arrays

      timeint_cdsp = 0

!EP   Initialize to 0

      do k=1,n3m
       q1_me(k)    =0.0d0
       q2_me(k)    =0.0d0
       q3_me(k)    =0.0d0
       dens_me(k)  =0.0d0
       q1_rms(k)   =0.0d0
       q2_rms(k)   =0.0d0
       q3_rms(k)   =0.0d0
       dens_rms(k) =0.0d0
       disste(k) = 0.0d0
       dissth(k) = 0.0d0
       densq3_me(k) = 0.0d0
      enddo

#endif
      return
      end
