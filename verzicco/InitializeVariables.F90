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
      use stat_arrays
      use decomp_2d
      use AuxiliaryRoutines
      implicit none
      
!-------------------------------------------------
! Arrays for grid making
!-------------------------------------------------

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

!-------------------------------------------------
! Arrays for density boundary conditions    
!-------------------------------------------------

      call AllocateReal2DArray(denbs,1,n2,1,n1)
      call AllocateReal2DArray(denbn,1,n2,1,n1)

!-------------------------------------------------
! Arrays for statistics    
!-------------------------------------------------

#ifdef STATS
      call AllocateReal1DArray(q1_me,1,n3m)
      call AllocateReal1DArray(q2_me,1,n3m)
      call AllocateReal1DArray(q3_me,1,n3m)

      call AllocateReal1DArray(q1_rms,1,n3m)
      call AllocateReal1DArray(q2_rms,1,n3m)
      call AllocateReal1DArray(q3_rms,1,n3m)

      call AllocateReal1DArray(dens_me,1,n3m)
      call AllocateReal1DArray(dens_rms,1,n3m)
      call AllocateReal1DArray(densq3_me,1,n3m)

#ifdef BALANCE
      call AllocateReal1DArray(disste,1,n3m)
      call AllocateReal1DArray(dissth,1,n3m)
#endif

#endif

      !-------------------------------------------------
      ! Arrays with ghost cells
      !-------------------------------------------------
      call AllocateReal3DArray(q1,1,n3,xstart(2)-1,xend(2)+1,xstart(3)-1,xend(3)+1)
      call AllocateReal3DArray(q2,1,n3,xstart(2)-1,xend(2)+1,xstart(3)-1,xend(3)+1)
      call AllocateReal3DArray(q3,1,n3,xstart(2)-1,xend(2)+1,xstart(3)-1,xend(3)+1)
      call AllocateReal3DArray(pr,1,n3,xstart(2)-1,xend(2)+1,xstart(3)-1,xend(3)+1)
      call AllocateReal3DArray(dens,1,n3,xstart(2)-1,xend(2)+1,xstart(3)-1,xend(3)+1)
      call AllocateReal3DArray(dphhalo,1,n3m,xstart(2)-1,xend(2)+1,xstart(3)-1,xend(3)+1)

      !-----------------------------------------------
      ! Arrays without ghost cells
      !-----------------------------------------------
      call AllocateReal3DArray(rhs,1,n3,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(dph,1,n3m,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(dq,1,n3,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(qcap,1,n3,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(ru1,1,n3,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(ru2,1,n3,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(ru3,1,n3,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(hro,1,n3,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(ruro,1,n3,xstart(2),xend(2),xstart(3),xend(3))

      return 
      end   

