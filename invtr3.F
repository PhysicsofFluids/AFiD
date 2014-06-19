c***********************************************************************
      subroutine invtr3
      use param
      use local_arrays, only: q3,rhs,ru3,qcap,pr
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,kc
      integer :: km,kp,ic
      real    :: alre,udx3
      real    :: amm,acc,app,dpx33,dq33

      alre=al/ren

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(xstart,xend,n3m,q3,pr)
!$OMP$ SHARED(kmv,kpv,am3ck,ac3ck,ap3ck)
!$OMP$ SHARED(al,ga,ro,alre,dt,qcap)
!$OMP$ SHARED(udx3c,rhs,ru3)
!$OMP$ PRIVATE(ic,jc,kc,km,kp)
!$OMP$ PRIVATE(amm,acc,app,udx3)
!$OMP$ PRIVATE(dq33,dpx33)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      do kc=2,n3m
      km=kc-1
      kp=kc+1
      udx3 = al*udx3c(kc)
      amm=am3ck(kc)
      acc=ac3ck(kc)
      app=ap3ck(kc)

c   33 second derivatives of q3
c
            dq33=q3(kp,jc,ic)*app
     %          +q3(kc,jc,ic)*acc
     %          +q3(km,jc,ic)*amm

c  component of grad(pr) along x3 direction
c
            dpx33=(pr(kc,jc,ic)-pr(km,jc,ic))*udx3
cm=======================================================     
            rhs(kc,jc,ic)=(ga*qcap(kc,jc,ic)+ro*ru3(kc,jc,ic)
     %                    +alre*dq33-dpx33)*dt 
cm=======================================================
c
c  updating of the non-linear terms
c
            ru3(kc,jc,ic)=qcap(kc,jc,ic)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

#ifdef SERIAL_DEBUG
      cksum=0.0d0
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      cksum=cksum+rhs(kc,jc,ic)
      end do
      end do
      end do
      write(*,*) 'hdnl3, rhs cksum', cksum

      do ic=1,n1m
      cksum=0.0d0
      do jc=1,n2m
      do kc=1,n3m
      cksum=cksum+ru3(kc,jc,ic)
      end do
      end do
      end do
      write(*,*) 'hdnl3, ru3 cksum', cksum
#endif

      call solq3k

      return
      end
c
