!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitializeVariables.F90                        !
!    CONTAINS: subroutine InitializeVariables             !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets to zero all    !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine InitializeVariables
      use param
      use local_arrays
      use decomp_2d, only: xstart,xend
      use AuxiliaryRoutines
      implicit none
      
      pr=0.d0
      rhs=0.d0
      ru1=0.d0
      ru2=0.d0
      ru3=0.d0
      ruro=0.d0
      dph=0.d0
      dphhalo=0.d0
      q1=0.d0
      q2=0.d0
      q3=0.d0
      dens=1.d0

      write(*,*) 'llego aca'
      call AllocateReal1DArray(tc,1,n1)
      call AllocateReal1DArray(tm,1,n1)
      call AllocateReal1DArray(ak1,1,n1)
      call AllocateReal1DArray(ao,1,n1)

      call AllocateReal1DArray(rc,1,n2)
      call AllocateReal1DArray(rm,1,n2)
      call AllocateReal1DArray(ak2,1,n2)
      call AllocateReal1DArray(ap,1,n2)

      call AllocateReal1DArray(zz,1,n3)
      call AllocateReal1DArray(zm,1,n3)
      call AllocateReal1DArray(g3rc,1,n3)
      call AllocateReal1DArray(g3rm,1,n3)

      call AllocateReal1DArray(udx3c,1,n3)
      call AllocateReal1DArray(udx3m,1,n3)

      call AllocateReal1DArray(ap3ck,1,n3)
      call AllocateReal1DArray(ac3ck,1,n3)
      call AllocateReal1DArray(am3ck,1,n3)

      call AllocateReal1DArray(ap3sk,1,n3)
      call AllocateReal1DArray(ac3sk,1,n3)
      call AllocateReal1DArray(am3sk,1,n3)

      call AllocateReal1DArray(ap3ssk,1,n3)
      call AllocateReal1DArray(ac3ssk,1,n3)
      call AllocateReal1DArray(am3ssk,1,n3)

      call AllocateReal1DArray(amphk,1,n3)
      call AllocateReal1DArray(acphk,1,n3)
      call AllocateReal1DArray(apphk,1,n3)
 
      call AllocateInt1dArray(kmc,1,n3)
      call AllocateInt1dArray(kpc,1,n3)
      call AllocateInt1dArray(kmv,1,n3)
      call AllocateInt1dArray(kpv,1,n3)

      call AllocateReal2DArray(denbs,1,n2,1,n1)
      call AllocateReal2DArray(denbn,1,n2,1,n1)

      return 
      end   

