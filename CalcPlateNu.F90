!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcPlateNu.F90                                !
!    CONTAINS: subroutine CalcPlateNu                     !
!                                                         ! 
!    PURPOSE: Calculate the Nusselt number at the top     !
!     and bottom plates and output to a file.             !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine CalcPlateNu
      use param
      use local_arrays, only: temp
      use mpih
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: j,i
      real ::  nuslow, nusupp
      real :: del,deln
  

      nuslow = 0.d0
      nusupp = 0.d0
      del  = 1.0/(xc(2)-xc(1))
      deln = 1.0/(xc(nx)-xc(nxm))

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,temp,del,deln) &
!$OMP   SHARED(nxm,nx) &
!$OMP   PRIVATE(i,j) &
!$OMP   REDUCTION(+:nuslow) &
!$OMP   REDUCTION(+:nusupp)
      do i=xstart(3),xend(3)
         do j=xstart(2),xend(2)
           nuslow = nuslow + (temp(1,j,i)-temp(2,j,i))*del
           nusupp = nusupp + (temp(nxm,j,i)-temp(nx,j,i))*deln
        enddo
      end do
!$OMP END PARALLEL DO

      nuslow = nuslow / (nzm*nym)
      nusupp = nusupp / (nzm*nym)

      call MpiSumRealScalar(nuslow)
      call MpiSumRealScalar(nusupp)

      if(ismaster) then
       open(97,file="nu_plate.out",status='unknown', &
        access='sequential',position='append')
       write(97,546) time, nuslow, nusupp
 546   format(4(1x,e14.6))
       close(97)
      endif

      return         
      end                                                               
