      subroutine CalcDissipationNu
      use mpih
      use param
      use local_arrays,only: q1,q2,q3,dens
      use decomp_2d, only: xstart,xend
      use stat_arrays

      implicit none
      integer :: i,j,k
      integer :: imm,ipp,jmm,jpp,kmm,kpp,kp
      real :: udx3_m,udx3_c
      real :: h11,h12,h13,h21,h22,h23,h31,h32,h33
      real :: nuth,nute,volt
      real :: udx1,udx2,dissipte,dissipth

      
      nute = 0.0d0
      nuth = 0.0d0

      udx1=dx1
      udx2=dx2
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Dissipation rates
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       	
!
!                                   1  |         | 2
!                   dissipation:  ---- | nabla  u|
!                                  Re  |         |
!

!$OMP  PARALLEL DO &
!$OMP  DEFAULT(none) &
!$OMP  SHARED(xstart,xend,g3rm,pra,q1,q2,q3,dens) &
!$OMP  SHARED(udx1,udx2,udx3m,udx3c) &
!$OMP  SHARED(nzm,nym,nxm,ren,pec) &
!$OMP  SHARED(kpv,kmv,zz,zm) &
!$OMP  SHARED(disste,dissth) &
!$OMP  PRIVATE(i,j,k,imm,ipp,jmm,jpp,kp,kpp,kmm) &
!$OMP  PRIVATE(udx3_m,udx3_c,dissipte,dissipth) &
!$OMP  PRIVATE(h11,h12,h13,h21,h22,h23,h31,h32,h33) &
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


       h11=(q1(k,j,ipp)-q1(k,j,i))*udx1
       h12=(q1(k,jpp,i)-q1(k,j,i))*udx2
       h13=(q1(kpp,j,i)-q1(k,j,i))*udx3m(k)

       h21=(q2(k,j,ipp)-q2(k,j,i))*udx1
       h22=(q2(k,jpp,i)-q2(k,j,i))*udx2
       h23=(q2(kpp,j,i)-q2(k,j,i))*udx3m(k)

       h31=(q3(k,j,ipp)-q3(k,j,i))*udx1
       h32=(q3(k,jpp,i)-q3(k,j,i))*udx2
       h33=(q3(kp,j,i)-q3(k,j,i))*udx3c(k)

       dissipte = 2.0*(h11**2+h22**2+h33**2)+ &
               (h21+h12)**2+(h31+h13)**2+(h32+h23)**2


       nute = nute + dissipte*g3rm(k)*pra

       h31=(dens(k,j,ipp)-dens(k,j,imm))*udx1*0.5
       h32=(dens(k,jpp,i)-dens(k,jmm,i))*udx2*0.5
       h33=(dens(kp,j,i)-dens(kmm,j,i))*udx3_c

       dissipth  = h31*h31 + h32*h32 + h33*h33 

       nuth = nuth+dissipth*g3rm(k)


!$OMP CRITICAL
       disste(k) =  disste(k) + dissipte / ren / real(nzm) / real(nym)
       dissth(k) =  dissth(k) + dissipth / pec / real(nzm) / real(nym)
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
