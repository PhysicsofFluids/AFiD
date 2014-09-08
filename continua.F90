!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: continua.F90                                   !
!    CONTAINS: subroutine continua                        !
!                                                         ! 
!    PURPOSE: Write down the full flow snapshot for       !
!     restarting the simulation at a later date           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine continua
      use param
      use local_arrays, only: q1,q2,q3,dens
      use decomp_2d, only: nrank
      character*30 :: filnam1

      filnam1 = trim('continua_dens.h5')
      call hdf_write(filnam1,dens)
      filnam1 = trim('continua_q1.h5')
      call hdf_write(filnam1,q1)
      filnam1 = trim('continua_q2.h5')
      call hdf_write(filnam1,q2)
      filnam1 = trim('continua_q3.h5')
      call hdf_write(filnam1,q3)
      
      if (nrank .eq. 0) then !EP only write once
       open(13,file='continua_grid.dat',status='unknown')
       rewind(13)                                                      
       write(13,*) n1,n2,n3
       write(13,*) rext,time
       write(13,*) istr3,str3
       close(13)
      endif

      end subroutine continua
