!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: divg.F90                                       !
!    CONTAINS: subroutine divg                            !
!                                                         ! 
!    PURPOSE: Compute the divergence of the velocity      !
!     at every point for the pressure correction          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalculateLocalDivergence
      use param
      use local_arrays, only: vz,vy,vx,dph
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,jp,kc,kp,ic,ip
      real    :: usdtal,dqcap   

      usdtal = 1.d0/(dt*al)

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,vz,vy,vx,dz,dy,udx3m,usdtal) &
!$OMP   SHARED(dph,nxm,xend) &
!$OMP   PRIVATE(ic,jc,kc,ip,jp,kp) &
!$OMP   PRIVATE(dqcap)
      do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
          jp=jc+1
            do kc=1,nxm
              kp=kc+1
              dqcap= (vz(kc,jc,ip)-vz(kc,jc,ic))*dz &
                    +(vy(kc,jp,ic)-vy(kc,jc,ic))*dy &
                    +(vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)
              dph(kc,jc,ic)=dqcap*usdtal
            enddo
         enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
