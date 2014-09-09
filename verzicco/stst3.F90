      subroutine stst3
      use param
      use local_arrays, only: dens,q1,q2,q3
      use stat3_param
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: i,j,m
      real,dimension(xstart(2):xend(2),xstart(3):xend(3)) :: &
     &      q3cc,q1cc,q2cc,denscc
      character*70 :: filnam
      character*1 :: charm

!EP   Slabs
!EP   cell center only q3

      do m=1,9
!$OMP  PARALLEL DO DEFAULT(SHARED) &
!$OMP   PRIVATE(i,j)
        do i=xstart(3),xend(3)
         do j=xstart(2),xend(2)
           q1cc(j,i) = q1(kslab(m),j,i)
           q2cc(j,i) = q2(kslab(m),j,i)
           q3cc(j,i) = (q3(kslab(m),j,i)+q3(kslab(m)+1,j,i))*0.5
           denscc(j,i) = dens(kslab(m),j,i)
          enddo
         enddo
!$OMP  END PARALLEL DO
      write(charm,28) m
   28 format(i1.1)
      filnam='slab'//charm//'q1_'
      call dump2dslab(q1cc,filnam)
      filnam='slab'//charm//'q2_'
      call dump2dslab(q2cc,filnam)
      filnam='slab'//charm//'q3_'
      call dump2dslab(q3cc,filnam)
      filnam='slab'//charm//'dens_'
      call dump2dslab(denscc,filnam)
      enddo

      return
      end subroutine stst3

      subroutine InitializeSlabDump
      use param
      use decomp_2d, only: nrank
      use stat3_param
      implicit none
      integer :: i,k,j
      real :: zmloc
      character(len=4) :: dummy

!EP   Read from stst3.in
      
      open(unit=19,file='stst3.in',status='old')
        read(19,301) dummy
        read(19,*) (zslab(i),i=2,9)
301     format(a4)                
      close(19)

!EP   Compute which kslab corresponds to which zslab
      
      kslab = 1
      
        do k=2,n3m
          zmloc=zm(k)
          do j=2,9
            if(zm(k).gt.zslab(j).and.zm(k-1).lt.zslab(j)) then
             kslab(j) = k
            endif
          enddo
        enddo


!EP   Write probe and slab locations
      
      if (nrank.eq.0) then
      open(unit=23,file='stst3locs.out',status='unknown')
        rewind(23)
        write(23,*) (kslab(i),i=1,9)
      close(23)
      endif

      return
      end subroutine InitializeSlabDump
      
      subroutine dump2dslab(var,filnam)
      USE param
      use mpih
      USE hdf5
      use decomp_2d, only: xstart,xend,nrank
      IMPLICIT none

      real, intent(in) :: var(xstart(2):xend(2) &
     &                  ,xstart(3):xend(3))

      integer hdf_error

      integer(HID_T) :: file_id
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: timespace
      integer(HID_T) :: memspace

      integer(HID_T) :: dset

      integer(HSIZE_T) :: dims(2)
      integer(HSIZE_T) :: dims2(1)

      integer(HID_T) :: plist_id
      integer(HSIZE_T), dimension(2) :: data_count  
      integer(HSSIZE_T), dimension(2) :: data_offset 

      integer :: comm, info
      integer :: ndims,ndims2

      real :: tprfi
      integer :: itime

      character*70,intent(in) :: filnam
      character*70 :: namfile
      character*8 :: ipfi

!RO   File writing part

      tprfi = 1/tpin
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i8.8)

      namfile=trim('./stst3/'//trim(filnam)//trim(ipfi)//'.h5')

!RO   Sort out MPI definitions

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Set offsets and element counts
   
      ndims = 2
      ndims2 = 1

      dims(1)=n2m
      dims(2)=n1m
      dims2(1)=1

      call h5screate_simple_f(ndims, dims,  &
     &                        filespace, hdf_error)

      data_count(1) = xend(2)-xstart(2)+1
      data_count(2) = xend(3)-xstart(3)+1

      data_offset(1) = xstart(2)-1
      data_offset(2) = xstart(3)-1

      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, &
     &    hdf_error)
 
      call h5pset_fapl_mpio_f(plist_id, comm, info, &
     &  hdf_error)

      call h5fcreate_f(namfile, H5F_ACC_TRUNC_F, file_id, &
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
     &   var(xstart(2):xend(2),xstart(3):xend(3)), dims,  &
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace,  &
     &   xfer_prp = plist_id)

      call h5pclose_f(plist_id, hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(memspace, hdf_error)
      call h5fclose_f(file_id, hdf_error)

!EP   Write time
      if(nrank.eq.0) then
      call h5fopen_f(namfile, H5F_ACC_RDWR_F, file_id, hdf_error)
      call h5screate_simple_f(ndims2, dims2, timespace, hdf_error) 
      call h5dcreate_f(file_id, 'time', H5T_NATIVE_DOUBLE, &
     &                timespace, dset, hdf_error)

      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, time, &
     &       dims2,hdf_error)

      call h5dclose_f(dset, hdf_error)

      call h5sclose_f(filespace, hdf_error)
      call h5fclose_f(file_id, hdf_error)
      endif

      return                                                          
      end subroutine dump2dslab
