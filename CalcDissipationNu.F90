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
      integer :: imm,ipp,jmm,jpp,kmm,kpp,kp
      real :: udx3_m,udx3_c
      real :: hxx,hxy,hxz,hyx,hyy,hyz,hzx,hzy,hzz
      real :: tx,ty,tz
      real :: nuth,nute,volt
      real :: udz,udy,dissipte,dissipth

      
      nute = 0.0d0
      nuth = 0.0d0

      udy=dy
      udz=dz
      
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(none) &
!$OMP  SHARED(xstart,xend,g3rm,pra,vz,vy,vx,temp) &
!$OMP  SHARED(udz,udy,udx3m,udx3c) &
!$OMP  SHARED(nzm,nym,nxm,ren,pec) &
!$OMP  SHARED(kpv,kmv,zz,zm) &
!$OMP  SHARED(disste,dissth) &
!$OMP  PRIVATE(i,j,k,imm,ipp,jmm,jpp,kp,kpp,kmm) &
!$OMP  PRIVATE(udx3_m,udx3_c,dissipte,dissipth) &
!$OMP  PRIVATE(hxx,hxy,hxz,hyx,hyy,hyz,hzx,hzy,hzz) &
!$OMP  PRIVATE(tx,ty,tz) &
!$OMP  REDUCTION(+:nute) &
!$OMP  REDUCTION(+:nuth)
      do i=xstart(3),xend(3)
       imm= i-1
       ipp= i+1
        do j=xstart(2),xend(2)
        jmm=j-1
        jpp=j+1

        do k=1,nxm
        kp=k+1
        kpp=kpv(k)
        kmm=kmv(k)

       udx3_m=1.0/(zm(kpp)-zm(kmm))
       udx3_c=1.0/(zz(kpp)-zz(kmm))

!
!      Viscous dissipation rate
!
!                       1  |         | 2
!                     ---- | nabla  u|
!                      Re  |         |
!

       hxx=(vx(kp,j,i)-vx(k,j,i))*udx3c(k)
       hxy=(vx(k,jpp,i)-vx(k,j,i))*udy
       hxz=(vx(k,j,ipp)-vx(k,j,i))*udz

       hyx=(vy(kpp,j,i)-vy(k,j,i))*udx3m(k)
       hyy=(vy(k,jpp,i)-vy(k,j,i))*udy
       hyz=(vy(k,j,ipp)-vy(k,j,i))*udz

       hzx=(vz(kpp,j,i)-vz(k,j,i))*udx3m(k)
       hzy=(vz(k,jpp,i)-vz(k,j,i))*udy
       hzz=(vz(k,j,ipp)-vz(k,j,i))*udz

       dissipte = 2.0*(hxx**2+hyy**2+hzz**2)+ &
               (hyz+hzy)**2+(hxz+hzx)**2+(hxy+hyx)**2


       nute = nute + dissipte*g3rm(k)*pra

!
!      Thermal gradient dissipation rate
!
!                       1  |         | 2
!                     ---- | nabla  T|
!                      Pe  |         |
!

       tx=(temp(kp,j,i )-temp(kmm,j,i))*udx3_c
       ty=(temp(k,jpp,i)-temp(k,jmm,i))*udy*0.5
       tz=(temp(k,j,ipp)-temp(k,j,imm))*udz*0.5

       dissipth  = tx*tx + ty*ty + tz*tz 

       nuth = nuth+dissipth*g3rm(k)


!$OMP CRITICAL
       disste(k) =  disste(k) + dissipte / (ren*real(nym)*real(nzm))
       dissth(k) =  dissth(k) + dissipth / (pec*real(nym)*real(nzm))
!$OMP END CRITICAL

       end do
       end do
       end do
!$OMP  END PARALLEL DO


       call MpiSumRealScalar(nuth)
       call MpiSumRealScalar(nute)
      
       volt = 1.d0/(real(nxm)*real(nzm)*real(nym))

      if(ismaster) then
      nute = nute*volt + 1
      nuth = nuth*volt 
      open(92,file='nu_diss.out',status='unknown',access='sequential', &
       position='append')
      write(92,*) time,nute,nuth
      close(92)
      endif

      return   
      end
