!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: divgck.F90                                     !
!    CONTAINS: subroutine divgck                          !
!                                                         ! 
!    PURPOSE: Check the maximum divergence of velocity    !
!     in the domain                                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine divgck(qmax)
      use param
      use local_arrays, only: q2,q3,q1
      use mpih
      use decomp_2d, only: xstart,xend,DECOMP_2D_COMM_CART_X
      implicit none
      real,intent(out) :: qmax
      integer :: jc,kc,kp,jp,ic,ip
      real    :: dqcap,my_qmax
        
!m==================
      my_qmax =-huge(0.0)
!m==================                                                                       
      qmax=0.d0                                                     
!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,n3m,q1,q2,q3,dx1,dx2,udx3m) &
!$OMP   PRIVATE(ic,jc,kc,ip,jp,kp) &
!$OMP   PRIVATE(dqcap) &
!$OMP   REDUCTION(max:my_qmax)
      do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
          jp=jc+1
            do kc=1,n3m
            kp=kc+1
              dqcap= (q1(kc,jc,ip)-q1(kc,jc,ic))*dx1 &
                    +(q2(kc,jp,ic)-q2(kc,jc,ic))*dx2 &
                    +(q3(kp,jc,ic)-q3(kc,jc,ic))*udx3m(kc)
              my_qmax = max(abs(dqcap),my_qmax)          
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
     
      call MPI_ALLREDUCE(my_qmax,qmax,1,MDP,MPI_MAX,DECOMP_2D_COMM_CART_X,ierr)
      
      return     
      end         
