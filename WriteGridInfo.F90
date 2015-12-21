!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: WriteGridInfo.F90                              !
!    CONTAINS: subroutine WriteGridInfo                   !
!                                                         ! 
!    PURPOSE: Write the grid information in               !
!     cordin_info.h5                                      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
       call HdfSerialWriteReal1D(dsetname,namfile,xm,nxm)
       dsetname = trim('xc')
       call HdfSerialWriteReal1D(dsetname,namfile,xc,nx)
       dsetname = trim('ym')
       call HdfSerialWriteReal1D(dsetname,namfile,ym,nym)
       dsetname = trim('zm')
       call HdfSerialWriteReal1D(dsetname,namfile,zm,nzm)

      endif

      return
      end


