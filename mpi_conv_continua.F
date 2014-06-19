      program mpi_conv_continua
      use hdf5
      implicit none
      integer m1m,m2m,m3,i,j,k
      parameter(m1m=32,m2m=32,m3=65)

      real, dimension(1:m1m,1:m2m,1:m3)::denso
      real, dimension(1:m1m,1:m2m,1:m3)::q1o
      real, dimension(1:m1m,1:m2m,1:m3)::q2o
      real, dimension(1:m1m,1:m2m,1:m3)::q3o

      real, dimension(1:m3,1:m2m,1:m1m)::dens
      real, dimension(1:m3,1:m2m,1:m1m)::q1
      real, dimension(1:m3,1:m2m,1:m1m)::q2
      real, dimension(1:m3,1:m2m,1:m1m)::q3

      real ckd, ckq1, ckq2, ckq3

      integer(HID_T) :: dset_q1
      integer(HID_T) :: dset_q2
      integer(HID_T) :: dset_q3
      integer(HID_T) :: dset_dens

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace

      integer(HID_T) :: dset_qua

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info
      integer :: ndims

      character*70 :: filnam1,filnam2,filnam3,filnam4



      call h5open_f(hdf_error)


      filnam1 = 'continua_dens.h5'
      filnam2 = 'continua_q1.h5'
      filnam3 = 'continua_q2.h5'
      filnam4 = 'continua_q3.h5'
   
      write(*,*) 'here'

!RO   Set offsets and element counts
      ndims = 3

      dims(1)=m1m
      dims(2)=m2m
      dims(3)=m3

      call h5fopen_f(filnam1, H5F_ACC_RDONLY_F, file_id,
     & hdf_error)

      call h5dopen_f(file_id, 'dens' ,dset_qua, hdf_error)

       call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE,
     &   denso(1:m1m,1:m2m,1:m3), dims, hdf_error)

      call h5dclose_f(dset_qua, hdf_error)

      call h5fclose_f(file_id, hdf_error)


!q1

      write(*,*) 'here'
      call h5fopen_f(filnam2, H5F_ACC_RDONLY_F, file_id,
     & hdf_error)

      call h5dopen_f(file_id, 'q1' ,dset_qua, hdf_error)

       call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE,
     &   q1o(1:m1m,1:m2m,1:m3), dims, hdf_error)

      call h5dclose_f(dset_qua, hdf_error)

      call h5fclose_f(file_id, hdf_error)

!q2

      write(*,*) 'here'
      call h5fopen_f(filnam3, H5F_ACC_RDONLY_F, file_id,
     & hdf_error)

      call h5dopen_f(file_id, 'q2' ,dset_qua, hdf_error)

       call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE,
     &   q2o(1:m1m,1:m2m,1:m3), dims, hdf_error)

      call h5dclose_f(dset_qua, hdf_error)

      call h5fclose_f(file_id, hdf_error)

!q3

      write(*,*) 'here'
      call h5fopen_f(filnam4, H5F_ACC_RDONLY_F, file_id,
     & hdf_error)

      call h5dopen_f(file_id, 'q3' ,dset_qua, hdf_error)

       call h5dread_f(dset_qua, H5T_NATIVE_DOUBLE,
     &   q3o(1:m1m,1:m2m,1:m3), dims, hdf_error)

      call h5dclose_f(dset_qua, hdf_error)

      call h5fclose_f(file_id, hdf_error)

!RO   Transpose

      do i=1,m1m
       do j=1,m2m
        do k=1,m3
         ckd=denso(i,j,k)+ckd
         ckq1=q1o(i,j,k)+ckq1
         ckq2=q2o(i,j,k)+ckq2
         ckq3=q3o(i,j,k)+ckq3
       end do
       end do
       end do

      write(*,*) ckd, ckq1, ckq2, ckq3
      write(*,*) 'transpose'

      do i=1,m1m
       do j=1,m2m
        do k=1,m3
         dens(k,j,i)=denso(i,j,k)
         q1(k,j,i)=q1o(i,j,k)
         q2(k,j,i)=q2o(i,j,k)
         q3(k,j,i)=q3o(i,j,k)
       end do
       end do
       end do

      ckd=0.0d0
      ckq1=0.0d0
      ckq2=0.0d0
      ckq3=0.0d0

      do i=1,m1m
       do j=1,m2m
        do k=1,m3
         ckd=dens(k,j,i)+ckd
         ckq1=q1(k,j,i)+ckq1
         ckq2=q2(k,j,i)+ckq2
         ckq3=q3(k,j,i)+ckq3
       end do
       end do
       end do

      write(*,*) ckd, ckq1, ckq2, ckq3

!RO   Set offsets and element counts
   
      dims(1)=m3
      dims(2)=m2m
      dims(3)=m1m

      call h5screate_simple_f(ndims, dims, filespace, hdf_error)

!EP   dens


      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_id,
     & hdf_error)


      call h5dcreate_f(file_id, 'dens', H5T_NATIVE_DOUBLE,
     &                filespace,dset_dens, hdf_error)


       call h5dwrite_f(dset_dens, H5T_NATIVE_DOUBLE,
     &   dens(1:m3,1:m2m,1:m1m), dims, hdf_error)

      call h5dclose_f(dset_dens, hdf_error)

      call h5fclose_f(file_id, hdf_error)

!EP   q1

      call h5fcreate_f(filnam2, H5F_ACC_TRUNC_F, file_id,
     & hdf_error)


      call h5dcreate_f(file_id, 'q1', H5T_NATIVE_DOUBLE,
     &                filespace, dset_q1, hdf_error)


       call h5dwrite_f(dset_q1, H5T_NATIVE_DOUBLE,
     &   q1(1:m3,1:m2m,1:m1m), dims, hdf_error)

      call h5dclose_f(dset_q1, hdf_error)
      call h5fclose_f(file_id, hdf_error)


!EP   q2

      call h5fcreate_f(filnam3, H5F_ACC_TRUNC_F, file_id,
     & hdf_error)


      call h5dcreate_f(file_id, 'q2', H5T_NATIVE_DOUBLE,
     &                filespace,dset_q2, hdf_error)


       call h5dwrite_f(dset_q2, H5T_NATIVE_DOUBLE,
     &   q2(1:m3,1:m2m,1:m1m), dims, hdf_error)

      call h5dclose_f(dset_q2, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   q3

      call h5fcreate_f(filnam4, H5F_ACC_TRUNC_F, file_id,
     & hdf_error)


      call h5dcreate_f(file_id, 'q3', H5T_NATIVE_DOUBLE,
     &                filespace, dset_q3, hdf_error)


       call h5dwrite_f(dset_q3, H5T_NATIVE_DOUBLE,
     &   q3(1:m3,1:m2m,1:m1m), dims, hdf_error)

      call h5dclose_f(dset_q3, hdf_error)
      call h5fclose_f(file_id, hdf_error)
      
      call h5sclose_f(filespace, hdf_error)

      end 

