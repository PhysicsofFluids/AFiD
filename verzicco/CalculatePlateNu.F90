
      subroutine CalculatePlateNu
      use param
      use local_arrays, only: dens
      use mpih
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: j,i
      real ::  anussupp,anusslow
      real :: del,deln
  

      anusslow = 0.d0
      anussupp = 0.d0
      del  = 1.0/(zz(2)-zz(1))
      deln = 1.0/(zz(nx)-zz(nxm))

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,dens,del,deln) &
!$OMP   SHARED(nxm,nx) &
!$OMP   PRIVATE(i,j) &
!$OMP   REDUCTION(+:anusslow) &
!$OMP   REDUCTION(+:anussupp)
      do i=xstart(3),xend(3)
         do j=xstart(2),xend(2)
           anusslow = anusslow + (dens(1,j,i)-dens(2,j,i))*del
           anussupp = anussupp + (dens(nxm,j,i)-dens(nx,j,i))*deln
        enddo
      end do
!$OMP END PARALLEL DO

      anusslow = anusslow / (nzm*nym)
      anussupp = anussupp / (nzm*nym)

      call MpiSumRealScalar(anusslow)
      call MpiSumRealScalar(anussupp)


      if(ismaster) then
       write(97,546) time, anusslow, anussupp
 546   format(4(1x,e14.6))
      endif

      return         
      end                                                               
