!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: hdf_write_serial_1d.F90                        !
!    CONTAINS: subroutine hdf_write_serial_1d             !
!                                                         ! 
!    PURPOSE: I/O routine. Write out a 1D array in        !
!     serial.                                             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hdf_write_serial_1d(dsetname,filename,n,var)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: n
      real, intent(in) :: var(n)
      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: dset
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)

!RO   Sort out MPI definitions and open file

      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdf_error)

!RO   Create dataspace

      dims(1)=n
      call h5screate_simple_f(1, dims, &
     &                        filespace, hdf_error)


!RO   Create the dataset with default properties.
      call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, &
     &                filespace, dset, hdf_error)

!RO   Select hyperslab  and then write it
       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
     &   var(1:n), dims, hdf_error)

!RO   Close properties and file

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(filespace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
      
      end subroutine hdf_write_serial_1d
