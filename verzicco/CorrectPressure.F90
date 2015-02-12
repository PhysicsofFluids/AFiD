!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: prcalc.F90                                     !
!    CONTAINS: subroutine prcalc                          !
!                                                         ! 
!    PURPOSE: Apply the pressure correction to the        !
!     pressure                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CorrectPressure
      use param
      use local_arrays, only: pr,dphhalo
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kp,km,jm,jp,jc,kc,ic,ip,im
      real    :: be,amm,acc,app

      be=al*beta
!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(pr,dphhalo,be,amphk,acphk,apphk) &
!$OMP   SHARED(xstart,xend,nxm,kmv,kpv,dzq,dyq) &
!$OMP   PRIVATE(ic,jc,kc) &
!$OMP   PRIVATE(im,jm,km,ip,jp,kp) &
!$OMP   PRIVATE(amm,acc,app)
      do ic=xstart(3),xend(3)
        im=ic-1
        ip=ic+1
        do jc=xstart(2),xend(2)
          jm=jc-1
          jp=jc+1
          do kc=1,nxm
            kp=kpv(kc)
            km=kmv(kc)
            amm=amphk(kc)
            acc=acphk(kc)
            app=apphk(kc)
              pr(kc,jc,ic)=pr(kc,jc,ic)+dphhalo(kc,jc,ic)-be*( &
              (dphhalo(kc,jc,ip) &
              -2.0*dphhalo(kc,jc,ic) &
              +dphhalo(kc,jc,im))*dzq+ &
              (dphhalo(kc,jp,ic) &
              -2.0*dphhalo(kc,jc,ic) &
              +dphhalo(kc,jm,ic))*dyq+ &
              (dphhalo(kp,jc,ic)*app &
              +dphhalo(kc,jc,ic)*acc &
              +dphhalo(km,jc,ic)*amm))
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
