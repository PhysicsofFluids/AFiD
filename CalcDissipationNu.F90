!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcDissipationNu.F90                          !
!    CONTAINS: subroutine CalcDissipationNu               !
!                                                         ! 
!    PURPOSE: Calculate the Nusselt number through the    !
!     global balance equations relating dissipation and   !
!     heat transport.                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcDissipationNu
      use mpih
      use param
      use local_arrays,only: vz,vy,vx,temp
      use decomp_2d, only: xstart,xend
      use stat_arrays

      implicit none
      integer :: i,j,k
      integer :: ip,jp,kp
      real :: dxm(nxm),dxc(nxm),lem(nxm),lec(nxm)
      real :: hxx,hxy,hxz,hyx,hyy,hyz,hzx,hzy,hzz
      real :: tx,ty,tz
      real :: nu_th,nu_mu,volt
      real :: dissipte,dissipte2,dissipth
      
      nu_mu = 0.0d0
      nu_th = 0.0d0
      
! RS: Calculated geometrical quantities outside the main loop.
! RS: Could potentially be moved to initialization routine 
      do k=1,nxm
        kp=k+1
        dxm(k)=1.0/(xm(kp)-xm(k))
        dxc(k)=1.0/(xc(kp)-xc(k))
        lem(k)=(xc(kp)-xc(k))*dx
        lec(k)=(xm(kp)-xm(k))*dx
      enddo
      lec(nxm)=(xm(nxm)-xm(nxm-1))*dx
      dxm(nxm)=1.0/(xm(nxm)-xm(nxm-1))

!$OMP  PARALLEL DO &
!$OMP  DEFAULT(none) &
!$OMP  SHARED(xstart,xend,vz,vy,vx,temp) &
!$OMP  SHARED(nxm,ren,pec,pra) &
!$OMP  SHARED(dxc,dxm,lec,lem) &
!$OMP  SHARED(disste,dissth) &
!$OMP  PRIVATE(i,j,k,ip,jp,kp) &
!$OMP  PRIVATE(dissipte,dissipte2,dissipth) &
!$OMP  PRIVATE(hxx,hxy,hxz,hyx,hyy,hyz,hzx,hzy,hzz) &
!$OMP  PRIVATE(tx,ty,tz) &
!$OMP  REDUCTION(+:nu_mu) &
!$OMP  REDUCTION(+:nu_th)
      do i=xstart(3),xend(3)
       ip= i+1
        do j=xstart(2),xend(2)
        jp=j+1

        do k=1,nxm
        kp=k+1

!       Viscous dissipation rate
!                       1  |         | 2
!                     ---- | nabla  u|
!                      Re  |         |
       hxx=(vx(kp,j,i)-vx(k,j,i))*dxc(k)
       hxy=(vx(k,jp,i)-vx(k,j,i))*dy
       hxz=(vx(k,j,ip)-vx(k,j,i))*dz

       hyx=(vy(kp,j,i)-vy(k,j,i))*dxm(k)
       hyy=(vy(k,jp,i)-vy(k,j,i))*dy
       hyz=(vy(k,j,ip)-vy(k,j,i))*dz

       hzx=(vz(kp,j,i)-vz(k,j,i))*dxm(k)
       hzy=(vz(k,jp,i)-vz(k,j,i))*dy
       hzz=(vz(k,j,ip)-vz(k,j,i))*dz

       dissipte  = (hxx*hxx+hxy*hxy+hxz*hxz) 
       dissipte2 = (hyx*hyx+hyy*hyy+hyz*hyz)+(hzx*hzx+hzy*hzy+hzz*hzz)

       nu_mu = nu_mu + dissipte*lem(k)*pra+dissipte2*lec(k)*pra 

!      Thermal gradient dissipation rate
!                       1  |         | 2
!                     ---- | nabla  T|
!                      Pe  |         |

       tx=(temp(kp,j,i)-temp(k,j,i))*dxc(k)
       ty=(temp(k,jp,i)-temp(k,j,i))*dy
       tz=(temp(k,j,ip)-temp(k,j,i))*dz

       dissipth  = tx*tx + ty*ty + tz*tz 
       nu_th = nu_th+dissipth*lem(k)

!$OMP CRITICAL
       disste(k) =  disste(k) + (dissipte + dissipte2) / (ren*real(nym)*real(nzm))
       dissth(k) =  dissth(k) + dissipth               / (pec*real(nym)*real(nzm))
!$OMP END CRITICAL

       end do
       end do
       end do
!$OMP  END PARALLEL DO

       call MpiSumRealScalar(nu_th)
       call MpiSumRealScalar(nu_mu)
      
       volt = 1.d0/(real(nxm)*real(nzm)*real(nym))

       if(ismaster) then
       nu_mu = nu_mu*volt + 1
       nu_th = nu_th*volt 
       open(92,file='nu_diss.out',status='unknown',access='sequential',position='append')
       write(92,*) time,nu_mu,nu_th
       close(92)
       endif

       return   
       end
