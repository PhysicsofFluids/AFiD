      subroutine invtr1
      use param
      use local_arrays, only: q1,dq,ru1,rhs,pr
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: kc,jc,ic,imm
      integer :: kmm,kpp
      real    :: alre,amm,acc,app
      real    :: d33q1,dpx11

c
      alre=al/ren


!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(xstart,xend,n3m,q1,pr)
!$OMP$ SHARED(kmv,kpv,am3sk,ac3sk,ap3sk)
!$OMP$ SHARED(dx1,al,ga,ro,alre,dt,dq)
!$OMP$ SHARED(udx3m,rhs,ru1)
!$OMP$ PRIVATE(ic,jc,kc,imm,kmm,kpp)
!$OMP$ PRIVATE(amm,acc,app)
!$OMP$ PRIVATE(d33q1,dpx11)
      do ic=xstart(3),xend(3)
      imm=ic-1
      do jc=xstart(2),xend(2)
      do kc=1,n3m
      kmm=kmv(kc)
      kpp=kpv(kc)
      amm=am3sk(kc)
      acc=ac3sk(kc)
      app=ap3sk(kc)

c
c   33 second derivative of q1
c
            d33q1=q1(kpp,jc,ic)*app
     %           +q1(kc,jc,ic)*acc
     %           +q1(kmm,jc,ic)*amm
      
c
c   component of grad(pr) along 2 direction
c
            dpx11=(pr(kc,jc,ic)-pr(kc,jc,imm))*dx1*al
c
c
c
            rhs(kc,jc,ic)=(ga*dq(kc,jc,ic)+ro*ru1(kc,jc,ic)
     %                    +alre*d33q1-dpx11)*dt

       

cm===========================================================
c
            ru1(kc,jc,ic)=dq(kc,jc,ic)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO


      call solq12k(q1(1:n3,xstart(2):xend(2),xstart(3):xend(3)),
     %     rhs(1:n3,xstart(2):xend(2),xstart(3):xend(3)))


      return
      end
c
