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
      use stat_arrays
      use AuxiliaryRoutines
      implicit none
      
      call DestroyReal1DArray(zc)
      call DestroyReal1DArray(zm)
      call DestroyReal1DArray(ak1)
      call DestroyReal1DArray(ao)

      call DestroyReal1DArray(yc)
      call DestroyReal1DArray(ym)
      call DestroyReal1DArray(ak2)
      call DestroyReal1DArray(ap)

      call DestroyReal1DArray(xc)
      call DestroyReal1DArray(xm)
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

      call DestroyReal2DArray(tempbp)
      call DestroyReal2DArray(temptp)

      call DestroyReal1DArray(vx_me)
      call DestroyReal1DArray(vy_me)
      call DestroyReal1DArray(vz_me)

      call DestroyReal1DArray(vx_rms)
      call DestroyReal1DArray(vy_rms)
      call DestroyReal1DArray(vz_rms)

      call DestroyReal1DArray(temp_me)
      call DestroyReal1DArray(temp_rms)
      call DestroyReal1DArray(tempvx_me)

      call DestroyReal3DArray(vx)
      call DestroyReal3DArray(vy)
      call DestroyReal3DArray(vz)
      call DestroyReal3DArray(temp)

      call DestroyReal3DArray(pr)
      call DestroyReal3DArray(rhs)

      call DestroyReal3DArray(dph)
      call DestroyReal3DArray(dphhalo)

      call DestroyReal3DArray(rux)
      call DestroyReal3DArray(ruy)
      call DestroyReal3DArray(ruz)
      call DestroyReal3DArray(rutemp)

      return 
      end   

