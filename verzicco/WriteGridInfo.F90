      subroutine WriteGridInfo
      use mpih
      use param
      use hdf5
      use decomp_2d, only: nrank

      IMPLICIT none

      integer hdf_error

      integer(HID_T) :: file_grid
      integer(HID_T) :: dset_grid
      integer(HID_T) :: dspace_grid

      integer(HSIZE_T) :: dims_grid(1)

      character*70 namfile


      if (nrank.eq.0) then
      namfile='cordin_info.h5'

!RO   Write the grid information normally

      call h5fcreate_f(namfile,H5F_ACC_TRUNC_F,file_grid,hdf_error)

      dims_grid(1)=n1m
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'tm', H5T_NATIVE_DOUBLE,dspace_grid, &
          dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, tm(1:n1m), &
              dims_grid, hdf_error)


      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n2m
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)

      call h5dcreate_f(file_grid, 'rm', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, rm(1:n2m), &
             dims_grid,hdf_error)

      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=n3m
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'zm', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, zm(1:n3m), &
              dims_grid, hdf_error)


      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)


      dims_grid(1)=1
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'Ray', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, ray, &
              dims_grid, hdf_error)


      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      dims_grid(1)=1
      call h5screate_simple_f(1, dims_grid, dspace_grid, hdf_error)
      call h5dcreate_f(file_grid, 'Pra', H5T_NATIVE_DOUBLE, &
                      dspace_grid, dset_grid, hdf_error)

      call h5dwrite_f(dset_grid, H5T_NATIVE_DOUBLE, pra,  &
              dims_grid, hdf_error)


      call h5dclose_f(dset_grid, hdf_error)
      call h5sclose_f(dspace_grid, hdf_error)

      call h5fclose_f(file_grid, hdf_error)
      endif

      return
      end


