      subroutine initia
      use param
      use local_arrays
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: j,k,i
      
cm=================================

!$OMP PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(xstart,xend,n3)
!$OMP$ SHARED(pr,dph,rhs,ru1,ru2,ru3,ruro)
!$OMP$ PRIVATE(i,j,k)
      do i=xstart(3),xend(3)
      do j=xstart(2),xend(2)
      do k=1,n3
      pr(k,j,i)=0.d0
      dph(k,j,i)=0.d0
      rhs(k,j,i)=0.d0
      ru1(k,j,i)=0.d0
      ru2(k,j,i)=0.d0
      ru3(k,j,i)=0.d0
      ruro(k,j,i)=0.d0
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
      
!$OMP PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(xstart,xend,n3,dens)
!$OMP$ SHARED(q1,q2,q3,n1m,n2m)
!$OMP$ PRIVATE(i,j,k)
      do i=max(xstart(3)-1,1),min(xend(3)+1,n1m)
      do j=max(xstart(2)-1,1),min(xend(2)+1,n2m)
      do k=1,n3
      q1(k,j,i)=0.d0
      q2(k,j,i)=0.d0
      q3(k,j,i)=0.d0
      dens(k,j,i)=1.d0
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
cm==================================         

      do i=1,n1
        tc(i) = 0.d0
        tm(i) = 0.d0
      end do
      do j=1,n2
        rc(j) = 0.d0
        rm(j) = 0.d0
      end do
      do k=1,n3
        am3ssk(k) = 0.d0
        ac3ssk(k) = 0.d0
        ap3ssk(k) = 0.d0
        zz(k) = 0.d0
        zm(k) = 0.d0
        g3rc(k) = 0.d0
        g3rm(k) = 0.d0
      end do
      return 
      end   
