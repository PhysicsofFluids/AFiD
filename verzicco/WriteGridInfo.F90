      subroutine WriteGridInfo
      use mpih
      use param
      use hdf5

      IMPLICIT none

      character*70 namfile
      character*30 :: dsetname


      if (ismaster) then 
       namfile='cordin_info.h5'
       call HdfCreateBlankFile(namfile)

       dsetname = trim('xm')
       call HdfSerialWriteReal1D(dsetname,namfile,tm,nzm)
       dsetname = trim('ym')
       call HdfSerialWriteReal1D(dsetname,namfile,rm,nym)
       dsetname = trim('zm')
       call HdfSerialWriteReal1D(dsetname,namfile,zm,nxm)
       dsetname = trim('zz')
       call HdfSerialWriteReal1D(dsetname,namfile,zz,nx)

      endif

      return
      end


