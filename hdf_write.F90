!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: hdf_write.F90                                  !
!    CONTAINS: subroutine hdf_write                       !
!                                                         ! 
!    PURPOSE: I/O routine. Write out a 3D array in        !
!     parallel of size n1o*n2o*n3o.                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine hdf_write(filnam1,qua)
      use param
      use mpih
      use decomp_2d, only: xstart,xend,nrank
      use hdf5
      
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      real, intent(in), dimension(1:n3,xstart(2)-1:xend(2)+1, &
     & xstart(3)-1:xend(3)+1)::qua

      character*30,intent(in) :: filnam1

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n3
      dims(2)=n2m
      dims(3)=n1m

      call h5screate_simple_f(ndims, dims,  &
     &                        filespace, hdf_error)

      data_count(1) = n3
      data_count(2) = xend(2)-xstart(2)+1
      data_count(3) = xend(3)-xstart(3)+1

      data_offset(1) = 0
      data_offset(2) = xstart(2)-1
      data_offset(3) = xstart(3)-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
          hdf_error)

      call h5pset_fapl_mpio_f(plist_id, comm, info, &
     &  hdf_error)

      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id, &
     & hdf_error, access_prp=plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dcreate_f(file_id, 'var', H5T_NATIVE_DOUBLE, &
     &                filespace, &
     &                dset, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
     &                      data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
     &                        hdf_error)
       call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, &
     &   qua(1:n3,xstart(2):xend(2),xstart(3):xend(3)), dims,  &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
      
      end subroutine hdf_write

