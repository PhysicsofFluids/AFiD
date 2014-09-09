!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DeallocateVariables.F90                        !
!    CONTAINS: subroutine DeallocateVariables             !
!                                                         ! 
!    PURPOSE: Finalization routine. Deallocates all       !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine DeallocateVariables
      use param
      use local_arrays
      use decomp_2d, only: xstart,xend
      use AuxiliaryRoutines
      implicit none
      
      call DestroyReal1DArray(tc)
      call DestroyReal1DArray(tm)
      call DestroyReal1DArray(ak1)
      call DestroyReal1DArray(ao)

      call DestroyReal1DArray(rc)
      call DestroyReal1DArray(rm)
      call DestroyReal1DArray(ak2)
      call DestroyReal1DArray(ap)

      call DestroyReal1DArray(zz)
      call DestroyReal1DArray(zm)
      call DestroyReal1DArray(g3rc)
      call DestroyReal1DArray(g3rm)

      call DestroyReal1DArray(udx3c)
      call DestroyReal1DArray(udx3m)

      call DestroyReal1DArray(ap3ck)
      call DestroyReal1DArray(ac3ck)
      call DestroyReal1DArray(am3ck)

      call DestroyReal1DArray(ap3sk)
      call DestroyReal1DArray(ac3sk)
      call DestroyReal1DArray(am3sk)

      call DestroyReal1DArray(ap3ssk)
      call DestroyReal1DArray(ac3ssk)
      call DestroyReal1DArray(am3ssk)

      call DestroyReal1DArray(amphk)
      call DestroyReal1DArray(acphk)
      call DestroyReal1DArray(apphk)
 
      call DestroyInt1dArray(kmc)
      call DestroyInt1dArray(kpc)
      call DestroyInt1dArray(kmv)
      call DestroyInt1dArray(kpv)

      call DestroyReal2DArray(denbs)
      call DestroyReal2DArray(denbn)

      return 
      end   

