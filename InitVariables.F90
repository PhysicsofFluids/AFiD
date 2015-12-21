!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: InitVariables.F90                              !
!    CONTAINS: subroutine InitVariables                   !
!                                                         ! 
!    PURPOSE: Initialization routine. Sets to zero all    !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine InitVariables
      use param
      use local_arrays
      use stat_arrays
      use decomp_2d
      use AuxiliaryRoutines
      implicit none
      
!-------------------------------------------------
! Arrays for grid making
!-------------------------------------------------

      call AllocateReal1DArray(zc,1,nz)
      call AllocateReal1DArray(zm,1,nz)
      call AllocateReal1DArray(ak1,1,nz)
      call AllocateReal1DArray(ao,1,nz)

      call AllocateReal1DArray(yc,1,ny)
      call AllocateReal1DArray(ym,1,ny)
      call AllocateReal1DArray(ak2,1,ny)
      call AllocateReal1DArray(ap,1,ny)

      call AllocateReal1DArray(xc,1,nx)
      call AllocateReal1DArray(xm,1,nx)
      call AllocateReal1DArray(g3rc,1,nx)
      call AllocateReal1DArray(g3rm,1,nx)

      call AllocateReal1DArray(udx3c,1,nx)
      call AllocateReal1DArray(udx3m,1,nx)

      call AllocateReal1DArray(ap3ck,1,nx)
      call AllocateReal1DArray(ac3ck,1,nx)
      call AllocateReal1DArray(am3ck,1,nx)

      call AllocateReal1DArray(ap3sk,1,nx)
      call AllocateReal1DArray(ac3sk,1,nx)
      call AllocateReal1DArray(am3sk,1,nx)

      call AllocateReal1DArray(ap3ssk,1,nx)
      call AllocateReal1DArray(ac3ssk,1,nx)
      call AllocateReal1DArray(am3ssk,1,nx)

      call AllocateReal1DArray(amphk,1,nx)
      call AllocateReal1DArray(acphk,1,nx)
      call AllocateReal1DArray(apphk,1,nx)
 
      call AllocateInt1dArray(kmc,1,nx)
      call AllocateInt1dArray(kpc,1,nx)
      call AllocateInt1dArray(kmv,1,nx)
      call AllocateInt1dArray(kpv,1,nx)

!-------------------------------------------------
! Arrays for temperature boundary conditions    
!-------------------------------------------------

      call AllocateReal2DArray(tempbp,1,ny,1,nz)
      call AllocateReal2DArray(temptp,1,ny,1,nz)

!-------------------------------------------------
! Arrays for statistics    
!-------------------------------------------------

      if (statcal) then
       call AllocateReal1DArray(vx_me,1,nxm)
       call AllocateReal1DArray(vy_me,1,nxm)
       call AllocateReal1DArray(vz_me,1,nxm)
 
       call AllocateReal1DArray(vx_rms,1,nxm)
       call AllocateReal1DArray(vy_rms,1,nxm)
       call AllocateReal1DArray(vz_rms,1,nxm)
 
       call AllocateReal1DArray(temp_me,1,nxm)
       call AllocateReal1DArray(temp_rms,1,nxm)
       call AllocateReal1DArray(tempvx_me,1,nxm)

       if (disscal) then
        call AllocateReal1DArray(disste,1,nxm)
        call AllocateReal1DArray(dissth,1,nxm)
       end if
      end if

      !-------------------------------------------------
      ! Arrays with ghost cells
      !-------------------------------------------------
      call AllocateReal3DArray(vy,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(vz,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(vx,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(pr,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(temp,1,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
      call AllocateReal3DArray(dphhalo,1,nxm,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

      !-----------------------------------------------
      ! Arrays without ghost cells
      !-----------------------------------------------
      call AllocateReal3DArray(rhs,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(dph,1,nxm,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(dq,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(qcap,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(rux,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(ruy,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(ruz,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(hro,1,nx,xstart(2),xend(2),xstart(3),xend(3))
      call AllocateReal3DArray(rutemp,1,nx,xstart(2),xend(2),xstart(3),xend(3))

      return 
      end   

