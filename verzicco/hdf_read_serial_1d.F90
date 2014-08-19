      subroutine hdf_read_serial_1d(dsetname,filename,n,var)
      use hdf5
      implicit none
      character*30,intent(in) :: dsetname,filename
      integer, intent(in) :: n
      real :: var(n)
      integer(HID_T) :: file_id
      integer(HID_T) :: dset
      integer :: hdf_error
      integer(HSIZE_T) :: dims(1)

!RO   Sort out MPI definitions and open file

      call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, hdf_error)

!RO   Create dataspace

      dims(1)=n

!RO   Create the dataset with default properties.
      call h5dopen_f(file_id, dsetname, dset, hdf_error)

!RO   Select hyperslab  and then read it
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, var(1:n), dims,hdf_error)

!RO   Close properties and file

      call h5dclose_f(dset, hdf_error)

      call h5fclose_f(file_id, hdf_error)
      
      end subroutine hdf_read_serial_1d
