      subroutine globalquantities
!EP   This routine calculates the maximum velocities, volume averaged nusselt number and other quantities required for sanity checks
      use param
      use local_arrays, only: q2,q3,q1,dens
      use decomp_2d, only: xstart,xend,DECOMP_2D_COMM_CART_X,nrank
      use mpih
      implicit none
      real :: my_vmax2,my_vmax3,my_vmax1,my_densm,my_denmax
      real :: my_denmin
      integer :: jc,kc,kp,ic
      real :: anusin,my_anusin,vol,q3cen,fac2,denscen
      real :: my_q2_rms_vol,my_q1_rms_vol,my_q3_rms_vol
      real :: my_q1q2q3_rms_vol,q2_rms_vol,q1_rms_vol
      real :: q3_rms_vol,q1q2q3_rms_vol,rradpr

!EP   Initialize
      my_vmax1=-huge(0.0d0)
      my_vmax2=-huge(0.0d0)
      my_vmax3=-huge(0.0d0)
      my_denmax=-huge(0.0d0)
      my_denmin=huge(0.0d0)
      my_densm=0.0d0
      my_anusin=0.d0 
      anusin=0.d0    
      my_q1_rms_vol = 0.0d0
      my_q2_rms_vol = 0.0d0
      my_q3_rms_vol = 0.0d0
      my_q1q2q3_rms_vol = 0.0d0
      vmax = 0.0d0
      vol = 1.d0/(alx3*dx3*real(n1m)*real(n2m))


!EP   Loop over volume
      do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
          do kc=1,n3m
          kp = kc + 1
          fac2 = g3rc(kc)
          my_vmax1 = max(my_vmax1,abs(q1(kc,jc,ic)))
          my_vmax2 = max(my_vmax2,abs(q2(kc,jc,ic)))
          my_vmax3 = max(my_vmax3,abs(q3(kc,jc,ic)))
          my_denmax = max(my_denmax,dens(kc,jc,ic))
          my_denmin = min(my_denmin,dens(kc,jc,ic))
          q3cen = (q3(kc,jc,ic)+q3(kp,jc,ic))*0.5d0
          denscen = (dens(kc,jc,ic)+dens(kp,jc,ic))*0.5d0
          my_anusin=my_anusin+denscen*q3cen*fac2
          my_densm=my_densm+denscen*fac2
          my_q1_rms_vol = my_q1_rms_vol + fac2*q1(kc,jc,ic)**2
          my_q2_rms_vol = my_q2_rms_vol + fac2*q2(kc,jc,ic)**2
          my_q3_rms_vol = my_q3_rms_vol + fac2*q3(kc,jc,ic)**2
          my_q1q2q3_rms_vol = my_q1q2q3_rms_vol + fac2*
     &    (q1(kc,jc,ic)**2+q2(kc,jc,ic)**2+q3(kc,jc,ic)**2)  
          enddo
        enddo
      enddo

!EP   Reduce
      call MPI_REDUCE(my_densm,densm,1,MDP,MPI_SUM,0,
     & DECOMP_2D_COMM_CART_X,ierr)
      call MPI_REDUCE(my_anusin,anusin,1,MDP,MPI_SUM,0,
     & DECOMP_2D_COMM_CART_X,ierr)
      call MPI_ALLREDUCE(my_vmax1,vmax(1),1,MDP,MPI_MAX,
     & DECOMP_2D_COMM_CART_X,ierr)
      call MPI_ALLREDUCE(my_vmax2,vmax(2),1,MDP,MPI_MAX,
     & DECOMP_2D_COMM_CART_X,ierr)
      call MPI_ALLREDUCE(my_vmax3,vmax(3),1,MDP,MPI_MAX,
     & DECOMP_2D_COMM_CART_X,ierr)
      call MPI_ALLREDUCE(my_denmax,denmax,1,MDP,MPI_MAX,
     & DECOMP_2D_COMM_CART_X,ierr)
      call MPI_ALLREDUCE(my_denmin,denmin,1,MDP,MPI_MIN,
     & DECOMP_2D_COMM_CART_X,ierr)
      call MPI_REDUCE(my_q1_rms_vol,q1_rms_vol,1,MDP,MPI_SUM,0,
     & MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q2_rms_vol,q2_rms_vol,1,MDP,MPI_SUM,0,
     & MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q3_rms_vol,q3_rms_vol,1,MDP,MPI_SUM,0,
     & MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(my_q1q2q3_rms_vol,q1q2q3_rms_vol,1,MDP,
     & MPI_SUM,0,MPI_COMM_WORLD,ierr)

!EP   Write
      if(nrank.eq.0) then
!EP   nusse.out
      anusin=1.d0 + dsqrt(pra*ray)*anusin*vol
      write(95,*) time, anusin

!EP   rms_vel.out
      rradpr=dsqrt(ray/pra)
      q1_rms_vol=dsqrt(q1_rms_vol*vol)*rradpr
      q2_rms_vol=dsqrt(q2_rms_vol*vol)*rradpr
      q3_rms_vol=dsqrt(q3_rms_vol*vol)*rradpr
      q1q2q3_rms_vol=dsqrt(q1q2q3_rms_vol*vol)*rradpr
       write(94,*) time,q1_rms_vol,q2_rms_vol,q3_rms_vol,
     & q1q2q3_rms_vol
      endif

      return   
      end     
