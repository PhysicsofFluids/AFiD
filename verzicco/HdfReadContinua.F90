!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: HdfReadContinua.F90                            !
!    CONTAINS: subroutine HdfReadContinua                 !
!                                                         ! 
!    PURPOSE: I/O routine. Read in a 3D array in          !
!     parallel of size n1o*n2o*n3o.                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine HdfReadContinua(n1o,n2o,n3o,xs2,xe2,xs3,xe3,intvar,qua)
      use mpih
      use param
      use hdf5
      implicit none

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset_qua

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer, intent(in) :: intvar,n1o,n2o,n3o
      integer, intent(in) :: xs2,xe2,xs3,xe3
      integer :: ndims
      real, dimension(1:n3o,xs2-lvlhalo:xe2+lvlhalo, /
       xs3-lvlhalo:xe3+lvlhalo),intent(out)::qua

      character*70 :: filnam1

!EP   Select file and dataset based on intvar

      select case (intvar)
        case (1)
          filnam1 = trim('continua_vx.h5')
        case (2)
          filnam1 = trim('continua_vy.h5')
        case (3)
          filnam1 = trim('continua_vz.h5')
        case (4)
          filnam1 = trim('continua_temp.h5')
      end select

!RO   Set offsets and element counts
   
      ndims = 3

      dims(1)=n3o
      dims(2)=n2o-1
      dims(3)=n1o-1

      data_count(1) = n3o
      data_count(2) = xe2-xs2+1
      data_count(3) = xe3-xs3+1

      data_offset(1) = 0
      data_offset(2) = xs2-1
      data_offset(3) = xs3-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
          hdf_error)

      call h5pset_fapl_mpio_f(plist_id, mpi_comm_world, &
     &  mpi_info_null, hdf_error)

      call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id, &
     & hdf_error, access_prp=plist_id)

      call h5dopen_f(file_id, 'var', &
     &                dset_qua, hdf_error)

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

      call h5dget_space_f(dset_qua, slabspace, hdf_error)
      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F, &
     &                      data_offset, data_count, hdf_error)
      call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf_error) 
      call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, &
     &                        hdf_error)
       call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE, &
     &   qua(1:n3o,xs2:xe2,xs3:xe3), dims,  &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
     &   xfer_prp = plist_id)
      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset_qua, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

      end subroutine HdfReadContinua

