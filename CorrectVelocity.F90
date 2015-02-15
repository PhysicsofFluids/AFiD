!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectVelocity.F90                            !
!    CONTAINS: subroutine CorrectVelocity                 !
!                                                         ! 
!    PURPOSE: Update velocities with the pressure         !
!     correction to enforce incompresibility              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CorrectVelocity
      use param
      use local_arrays, only: vy,vx,dphhalo,vz
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,jm,kc,km,ic,im
      real    :: usukm,udy,udz,locdph

      udy = al*dt*dy
      udz = al*dt*dz

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(vz,vy,vx,dphhalo,udz,udy,udx3c) &
!$OMP   SHARED(xstart,xend,nxm,kmv,dt,al) &
!$OMP   PRIVATE(ic,jc,kc) &
!$OMP   PRIVATE(im,jm,km,usukm,locdph)
      do ic=xstart(3),xend(3)
        im=ic-1
        do jc=xstart(2),xend(2)
          jm=jc-1
          do kc=1,nxm
          km=kmv(kc)
          usukm = al*dt*udx3c(kc)
          locdph=dphhalo(kc,jc,ic)
          vx(kc,jc,ic)=vx(kc,jc,ic)- &
            (locdph-dphhalo(km,jc,ic))*usukm
          vy(kc,jc,ic)=vy(kc,jc,ic)- &
            (locdph-dphhalo(kc,jm,ic))*udy
          vz(kc,jc,ic)=vz(kc,jc,ic)- &
            (locdph-dphhalo(kc,jc,im))*udz
        enddo 
       enddo
      enddo
!$OMP END PARALLEL DO

      return
      end

