!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: divgck.F90                                     !
!    CONTAINS: subroutine divgck                          !
!                                                         ! 
!    PURPOSE: Check the maximum divergence of velocity    !
!     in the domain                                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CheckDivergence(qmax)
      use param
      use local_arrays, only: vy,vx,vz
      use mpih
      use decomp_2d, only: xstart,xend
      implicit none
      real,intent(out) :: qmax
      integer :: jc,kc,kp,jp,ic,ip
      real    :: dqcap
        
      qmax =-huge(0.0)

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vz,vy,vx,dz,dy,udx3m) &
!$OMP   PRIVATE(ic,jc,kc,ip,jp,kp) &
!$OMP   PRIVATE(dqcap) &
!$OMP   REDUCTION(max:qmax)
      do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
          jp=jc+1
            do kc=1,nxm
            kp=kc+1
              dqcap= (vz(kc,jc,ip)-vz(kc,jc,ic))*dz &
                    +(vy(kc,jp,ic)-vy(kc,jc,ic))*dy &
                    +(vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)
              qmax = max(abs(dqcap),qmax)          
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      call MpiMaxRealScalar(qmax)
     
      
      return     
      end         
