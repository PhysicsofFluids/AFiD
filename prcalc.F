
************************************************************************
c   this subroutine performs the calculation of the pressure.
c   this depends on the fractional step
c
      subroutine prcalc
      use param
      use local_arrays, only: pr,dph
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kp,km,jm,jp,jc,kc,ic,ip,im
      real    :: be,amm,acc,app
c
c    the pressure is evaluated at the center of the box.
c
c     p^{n+1} = p^{n} + phi^{n+1} - b * Nabla^2 phi^{n+1}
c
      be=al*beta
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(pr,dph,be,amphk,acphk,apphk)
!$OMP$ SHARED(xstart,xend,n3m,kmv,kpv,dx1q,dx2q)
!$OMP$ PRIVATE(ic,jc,kc)
!$OMP$ PRIVATE(im,jm,km,ip,jp,kp)
!$OMP$ PRIVATE(amm,acc,app)
      do ic=xstart(3),xend(3)
        im=ic-1
        ip=ic+1
        do jc=xstart(2),xend(2)
          jm=jc-1
          jp=jc+1
          do kc=1,n3m
            kp=kpv(kc)
            km=kmv(kc)
            amm=amphk(kc)
            acc=acphk(kc)
            app=apphk(kc)
              pr(kc,jc,ic)=pr(kc,jc,ic)+dph(kc,jc,ic)-be*(
     %        (dph(kc,jc,ip)
     %        -2.0*dph(kc,jc,ic)
     %        +dph(kc,jc,im))*dx1q+
     %        (dph(kc,jp,ic)
     %        -2.0*dph(kc,jc,ic)
     %        +dph(kc,jm,ic))*dx2q+
     %        (dph(kp,jc,ic)*app
     %        +dph(kc,jc,ic)*acc
     %        +dph(km,jc,ic)*amm))
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
      return
      end
c
c
