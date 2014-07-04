c***********************************************************************
      subroutine invtrro
      use param
      use local_arrays, only: dens,hro,ruro,rhs
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,kc,ic
      integer :: km,kp
      real    :: alpec,dq33
      real    :: app,acc,amm
c
c     h term for the q3 momentum equation at i+1/2,j+1/2,k
c

      alpec=al/pec

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(xstart,xend,n3m,dens)
!$OMP$ SHARED(kmv,kpv,am3ck,ac3ck,ap3ck)
!$OMP$ SHARED(ga,ro,alpec,dt)
!$OMP$ SHARED(rhs,ruro,hro)
!$OMP$ PRIVATE(ic,jc,kc,km,kp)
!$OMP$ PRIVATE(amm,acc,app)
!$OMP$ PRIVATE(dq33)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      do kc=2,n3m
      km=kc-1
      kp=kc+1
      app=ap3ck(kc)
      acc=ac3ck(kc)
      amm=am3ck(kc)

c   33 second derivatives of dens
c
               dq33= dens(kp,jc,ic)*app
     %              +dens(kc,jc,ic)*acc
     %              +dens(km,jc,ic)*amm
c     
c
c    right hand side of the density equation
c
cm===========================================================
            rhs(kc,jc,ic)=(ga*hro(kc,jc,ic)+ro*ruro(kc,jc,ic)
     %              +alpec*dq33)*dt
cm===========================================================
c
c    updating of the non-linear terms
c
            ruro(kc,jc,ic)=hro(kc,jc,ic)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO


      call solrok


!      write(*,*) dens(1,:,:)
       dens(1,xstart(2):xend(2),xstart(3):xend(3)) 
     4    = denbs(xstart(2):xend(2),xstart(3):xend(3))
       dens(n3,xstart(2):xend(2),xstart(3):xend(3)) 
     4    = denbn(xstart(2):xend(2),xstart(3):xend(3))
      
      return
      end
