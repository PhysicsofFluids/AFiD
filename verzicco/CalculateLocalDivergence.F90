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
      use local_arrays, only: q1,q2,q3,dph
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,jp,kc,kp,ic,ip
      real    :: usdtal,dqcap   

      usdtal = 1.d0/(dt*al)

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,q1,q2,q3,dx1,dx2,udx3m,usdtal) &
!$OMP   SHARED(dph,n3m,xend) &
!$OMP   PRIVATE(ic,jc,kc,ip,jp,kp) &
!$OMP   PRIVATE(dqcap)
      do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
          jp=jc+1
            do kc=1,n3m
              kp=kc+1
              dqcap= (q1(kc,jc,ip)-q1(kc,jc,ic))*dx1 &
                    +(q2(kc,jp,ic)-q2(kc,jc,ic))*dx2 &
                    +(q3(kp,jc,ic)-q3(kc,jc,ic))*udx3m(kc)
              dph(kc,jc,ic)=dqcap*usdtal
            enddo
         enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
