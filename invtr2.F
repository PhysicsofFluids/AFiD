c***********************************************************************
      subroutine invtr2
      use param
      use local_arrays, only: q2,ru2,pr,rhs,dph
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,jmm,jc,ic
      integer :: kpp,kmm
      real    :: alre,udx2
      real    :: amm,acc,app
      real    :: d33q2,dpx22


      alre=al/ren
      udx2=dx2*al

c
c     h term for the q2 momentum equation at i+1/2,j,k+1/2
c
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(xstart,xend,n3m,q2,pr)
!$OMP$ SHARED(kmv,kpv,am3sk,ac3sk,ap3sk)
!$OMP$ SHARED(dx2,al,ga,ro,alre,dt,dph)
!$OMP$ SHARED(udx2,udx3m,rhs,ru2)
!$OMP$ PRIVATE(ic,jc,kc,kmm,kpp,jmm)
!$OMP$ PRIVATE(amm,acc,app)
!$OMP$ PRIVATE(d33q2,dpx22)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      jmm=jc-1
      do kc=1,n3m
      kmm=kmv(kc)
      kpp=kpv(kc)
      amm=am3sk(kc)
      acc=ac3sk(kc)
      app=ap3sk(kc)

c
c   33 second derivative of q2
c
            d33q2=q2(kpp,jc,ic)*app
     %           +q2(kc,jc,ic)*acc
     %           +q2(kmm,jc,ic)*amm

c
c   component of grad(pr) along 2 direction
c
            dpx22=(pr(kc,jc,ic)-pr(kc,jmm,ic))*udx2

c
            rhs(kc,jc,ic)=(ga*dph(kc,jc,ic)+ro*ru2(kc,jc,ic)
     %                    +alre*d33q2-dpx22)*dt

cm===========================================================
c
            ru2(kc,jc,ic)=dph(kc,jc,ic)
      enddo
      enddo
      enddo

#ifdef SERIAL_DEBUG
      cksum=0.0d0
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      cksum=cksum+rhs(kc,jc,ic)
      end do
      end do
      end do
      write(*,*) 'hdnl2, rhs cksum', cksum

      cksum=0.0d0
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      cksum=cksum+ru2(kc,jc,ic)
      end do
      end do
      end do
      write(*,*) 'hdnl2, ru2 cksum', cksum
#endif

      call solq12k(q2(1:n3,xstart(2):xend(2),xstart(3):xend(3)),
     %     rhs(1:n3,xstart(2):xend(2),xstart(3):xend(3)))

#ifdef SERIAL_DEBUG
      cksum=0.0d0
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      cksum=cksum+q2(kc,jc,ic)
      end do
      end do
      end do
      write(*,*) 'hdnl2, q2 cksum', cksum
#endif
      
      return
      end
c
