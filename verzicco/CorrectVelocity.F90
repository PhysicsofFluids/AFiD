!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: updvp.F90                                      !
!    CONTAINS: subroutine updvp                           !
!                                                         ! 
!    PURPOSE: Update velocities with the pressure         !
!     correction to enforce incompresibility              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CorrectVelocity
      use param
      use local_arrays, only: q2,q3,dphhalo,q1
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,jm,kc,km,ic,im
      real    :: usukm,udx2,udx1,locdph

      udx1 = al*dt*dx1
      udx2 = al*dt*dx2

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(q1,q2,q3,dphhalo,udx1,udx2,udx3c) &
!$OMP   SHARED(xstart,xend,n3m,kmv,dt,al) &
!$OMP   PRIVATE(ic,jc,kc) &
!$OMP   PRIVATE(im,jm,km,usukm,locdph)
      do ic=xstart(3),xend(3)
        im=ic-1
        do jc=xstart(2),xend(2)
          jm=jc-1
          do kc=1,n3m
          km=kmv(kc)
          usukm = al*dt*udx3c(kc)
          locdph=dphhalo(kc,jc,ic)
          q1(kc,jc,ic)=q1(kc,jc,ic)- &
            (locdph-dphhalo(kc,jc,im))*udx1
          q2(kc,jc,ic)=q2(kc,jc,ic)- &
            (locdph-dphhalo(kc,jm,ic))*udx2
          q3(kc,jc,ic)=q3(kc,jc,ic)- &
            (locdph-dphhalo(km,jc,ic))*usukm
        enddo 
       enddo
      enddo
!$OMP END PARALLEL DO

      return
      end

