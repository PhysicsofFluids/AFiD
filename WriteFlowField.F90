!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: WriteFlowField.F90                             !
!    CONTAINS: subroutine WriteFlowField                  !
!                                                         ! 
!    PURPOSE: Write down the full flow snapshot for       !
!     restarting the simulation at a later date           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine WriteFlowField
      use param
      use local_arrays, only: q1,q2,q3,temp
      implicit none
      character*30 :: filnam1,dsetname

      filnam1 = trim('continua_temp.h5')
      call HdfWriteRealHalo3D(filnam1,temp)
      filnam1 = trim('continua_q1.h5')
      call HdfWriteRealHalo3D(filnam1,q1)
      filnam1 = trim('continua_q2.h5')
      call HdfWriteRealHalo3D(filnam1,q2)
      filnam1 = trim('continua_q3.h5')
      call HdfWriteRealHalo3D(filnam1,q3)

      if (ismaster) then !EP only write once
       filnam1 = trim('continua_master.h5')
       call HdfCreateBlankFile(filnam1)
 
       dsetname = trim('nx')
       call HdfSerialWriteIntScalar(dsetname,filnam1,nx)
       dsetname = trim('ny')
       call HdfSerialWriteIntScalar(dsetname,filnam1,ny)
       dsetname = trim('nz')
       call HdfSerialWriteIntScalar(dsetname,filnam1,nz)
       dsetname = trim('ylen')
       call HdfSerialWriteRealScalar(dsetname,filnam1,ylen)
       dsetname = trim('zlen')
       call HdfSerialWriteRealScalar(dsetname,filnam1,zlen)
       dsetname = trim('time')
       call HdfSerialWriteRealScalar(dsetname,filnam1,time)
       dsetname = trim('istr3')
       call HdfSerialWriteIntScalar(dsetname,filnam1,istr3)
       dsetname = trim('str3')
       call HdfSerialWriteRealScalar(dsetname,filnam1,str3)

      endif

      end subroutine WriteFlowField
