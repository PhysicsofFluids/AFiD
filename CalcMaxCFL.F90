
      subroutine CalcMaxCFL(cflm)
      use param
      use local_arrays, only: q2,q3,q1
      use decomp_2d
      use mpih
      implicit none
      real,intent(out)    :: cflm
      integer :: j,k,jp,kp,i,ip
      real :: qcf,udx3
      
      cflm=0.00000001d0
!                                                                       
!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,n3m,q1,q2,q3) &
!$OMP   SHARED(dx1,dx2,udx3m) &
!$OMP   PRIVATE(i,j,k,ip,jp,kp,udx3,qcf) &
!$OMP   REDUCTION(max:cflm)
      do i=xstart(3),xend(3)
        ip=i+1
        do j=xstart(2),xend(2)
          jp=j+1
          do k=1,n3m
           udx3=udx3m(k)
           kp=k+1
            qcf=( abs((q1(k,j,i)+q1(k,j,ip))*0.5d0*dx1) &
                 +abs((q2(k,j,i)+q2(k,jp,i))*0.5d0*dx2) &
                 +abs((q3(k,j,i)+q3(kp,j,i))*0.5d0*udx3))

            cflm = max(cflm,qcf)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
            
      call MpiMaxRealScalar(cflm)

      return  
      end
