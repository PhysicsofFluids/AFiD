      subroutine vmaxv
!EP   This routine calculates the maximum velocities and volume averaged nusselt number
      use param
      use local_arrays, only: q2,q3,q1,dens
      use decomp_2d, only: xstart,xend,DECOMP_2D_COMM_CART_X,nrank
      use mpih
      implicit none
      real    :: my_vmax2,my_vmax3,my_vmax1
      integer :: jc,kc,kp,ic
      real :: anusin,my_anusin,vol,q3cen,fac2,denscen

        my_vmax1=-100.d0
        my_vmax2=-100.d0
        my_vmax3=-100.d0
      my_anusin=0.d0 
      anusin=0.d0    
      vol = 1.d0/(alx3*dx3*real(n1m)*real(n2m))
      do ic=xstart(3),xend(3)
         do jc=xstart(2),xend(2)
            do kc=1,n3m
            kp = kc + 1
            fac2 = g3rc(kc)
              my_vmax1 = max(my_vmax1,abs(q1(kc,jc,ic)))
              my_vmax2 = max(my_vmax2,abs(q2(kc,jc,ic)))
              my_vmax3 = max(my_vmax3,abs(q3(kc,jc,ic)))
              q3cen = (q3(kc,jc,ic)+q3(kp,jc,ic))*0.5d0
              denscen = (dens(kc,jc,ic)+dens(kp,jc,ic))*0.5d0
              my_anusin=my_anusin+denscen*q3cen*fac2
       enddo
       enddo
       enddo

      call MPI_REDUCE(my_anusin,anusin,1,MDP,MPI_SUM,0,
     & DECOMP_2D_COMM_CART_X,ierr)
      call MPI_ALLREDUCE(my_vmax1,vmax(1),1,MDP,MPI_MAX,
     & DECOMP_2D_COMM_CART_X,ierr)
      call MPI_ALLREDUCE(my_vmax2,vmax(2),1,MDP,MPI_MAX,
     & DECOMP_2D_COMM_CART_X,ierr)
      call MPI_ALLREDUCE(my_vmax3,vmax(3),1,MDP,MPI_MAX,
     & DECOMP_2D_COMM_CART_X,ierr)

      if(nrank.eq.0) then
      anusin=1.d0 + dsqrt(pra*ray)*anusin*vol
      write(95,*) time, anusin
      endif

      return   
      end     
