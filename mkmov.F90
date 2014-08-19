!***********************************************************************
      subroutine mkmov
#ifdef MOVIE
      use local_arrays, only: dens
      use mpi_param
      use mpih
      use hdf5
      use param

      IMPLICIT NONE

      integer ic,jc,kc

      integer hdf_error
      integer(HID_T) :: filespace
      integer(HID_T) :: slabspace
      integer(HID_T) :: memspace

      integer(HID_T) :: file_densv

      integer(HID_T) :: dset_densv

      integer(HSIZE_T) :: dims(3)

      integer(HID_T) :: file_plist
      integer(HID_T) :: slab_plist
      integer(HSIZE_T), dimension(3) :: data_count  
      integer(HSSIZE_T), dimension(3) :: data_offset 

      integer :: comm, info

      integer i,j,k
      integer ndims,itime

      real :: tprfi
      character*70 filnam1,filnamxdm
      character*5 ipfi

!RO   Sort out MPI definitions and file names

      tprfi = 1/tframe
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i5.5)

      filnam1='movie/frame'//ipfi//'.h5'
      filnamxdm = './movie/frame'//ipfi//'.xmf' 

      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL

!RO   Set offsets and element counts

      ndims=3

      dims(1)=n1m
      dims(2)=n2m
      dims(3)=n3m

      data_count(1) = n1m
      data_count(2) = n2m
      data_count(3) = kend-kstart+1

      data_offset(1) = 0
      data_offset(2) = 0
      data_offset(3) = kstart-1 


!RO   Set up MPI file properties

      call h5pcreate_f(H5P_FILE_ACCESS_F, file_plist, hdf_error)
      call h5pset_fapl_mpio_f(file_plist, comm, info, hdf_error)

      call h5pcreate_f(H5P_DATASET_XFER_F, slab_plist, hdf_error) 
      call h5pset_dxpl_mpio_f(slab_plist, H5FD_MPIO_COLLECTIVE_F,
     &                        hdf_error)

!RO   Create dataspace

      call h5screate_simple_f(ndims, dims, 
     &                        filespace, hdf_error)

!RO   Create dataspace in memory

      call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 

!RO   Open first continua file for dens

      call h5fcreate_f(filnam1, H5F_ACC_TRUNC_F, file_densv, hdf_error,
     &                 access_prp=file_plist)

!RO   Create dataset on file

      call h5dcreate_f(file_densv, 'dens', H5T_NATIVE_DOUBLE,
     &                filespace, dset_densv, hdf_error)

!RO   Set hyperslab

      call h5dget_space_f(dset_densv, slabspace, hdf_error)

      call h5sselect_hyperslab_f (slabspace, H5S_SELECT_SET_F,
     &                      data_offset, data_count, hdf_error)

      call h5dwrite_f(dset_densv, H5T_NATIVE_DOUBLE,
     &   dens(1:n1m,1:n2m,kstart:kend), dims, 
     &   hdf_error, file_space_id = slabspace, mem_space_id = memspace, 
     &   xfer_prp = slab_plist)

!RO   Close dataset and file for dens

      call h5dclose_f(dset_densv, hdf_error)
      call h5fclose_f(file_densv, hdf_error)

!RO   Close all other stuff

      call h5sclose_f(memspace, hdf_error)
      call h5sclose_f(slabspace, hdf_error)
      call h5sclose_f(filespace, hdf_error)
      call h5pclose_f(file_plist, hdf_error)
      call h5pclose_f(slab_plist, hdf_error)

!EP   Write the xdm

      if (myid.eq.0) then

      open(45,file=filnamxdm,status='unknown')
      rewind(45)
      write(45,'("<?xml version=""1.0"" ?>")')
      write(45,'("<!DOCTYPE Xdmf SYSTEM ""Xdmf.dtd"" []>")')
      write(45,'("<Xdmf Version=""2.0"">")')
      write(45,'("<Domain>")')
      write(45,'("<Grid Name=""RB Cartesian"" GridType=""Uniform"">")')
      write(45,'("<Topology TopologyType=""3DRectMesh"" 
     &NumberOfElements=""",i4," ",i4," ",i4,"""/>")') n3m,n2m,n1m
      write(45,'("<Geometry GeometryType=""VXVYVZ"">")')
      write(45,'("<DataItem Dimensions=""",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n1m
      write(45,'("cordin_info.h5:/tm")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n2m
      write(45,'("cordin_info.h5:/rm")')
      write(45,'("</DataItem>")')
      write(45,'("<DataItem Dimensions=""",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')n3m
      write(45,'("cordin_info.h5:/zm")')
      write(45,'("</DataItem>")')
      write(45,'("</Geometry>")')
      write(45,'("<Attribute Name=""Temperature""
     & AttributeType=""Scalar"" Center=""Node"">")')
      write(45,'("<DataItem Dimensions=""",i4," ",i4," ",i4,"""
     & NumberType=""Float"" Precision=""4"" Format=""HDF"">")')
     &n3m,n2m,n1m
      write(45,'("frame",i5.5,".h5:/dens")') itime
      write(45,'("</DataItem>")')
      write(45,'("</Attribute>")')
      write(45,'("<Time Value=""",e12.5""" />")') time
      write(45,'("</Grid>")')
      write(45,'("</Domain>")')
      write(45,'("</Xdmf>")')
      close(45)

      endif

#endif
      return                                                          
      end                                                             

