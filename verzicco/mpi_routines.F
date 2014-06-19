c===============================================
      subroutine block(n, p, irank, istart, iend, blcsz)
      implicit none
      integer,intent(in) :: n,p,irank
      integer,intent(out) :: istart,iend
      integer :: i
      integer,dimension(0:p-1),intent(out) :: blcsz
      
      do i=0,p-1
      blcsz(i) = floor(real((n+p-i-1)/p))
      enddo
      istart = sum(blcsz(0:irank))-blcsz(irank)+1
      iend = istart+blcsz(irank)-1

      end subroutine block
c=================================================           
      subroutine mpi_workdistribution
      use param
      use mpih 
      use mpi_param
      use decomp_2d
      implicit none
      

!EP   Initialize pencil decomposition with automatic distribution and periodic sidewalls

      !------------------------------------------------
      
      end subroutine mpi_workdistribution 
      !------------------------------------------------
      !------------------------------------------------
      !------------------------------------------------
      subroutine mem_alloc
      use param
      use mpi_param
      use local_arrays
      use decomp_2d
#ifdef STATS
      use stat_arrays
#endif
      implicit none
!     type(DECOMP_INFO) :: decomp_halo
!     
!     decomp_halo = ph
!     
!EP   Make halo
!     decomp_halo%xst(1) = ph%xst(1) 
!     decomp_halo%xen(1) = ph%xen(1) 
!     decomp_halo%xst(2) = ph%xst(2) - 1
!     decomp_halo%xen(2) = ph%xen(2) + 1
!     decomp_halo%xst(3) = ph%xst(3) - 1
!     decomp_halo%xen(3) = ph%xen(3) + 1
      

      !-------------------------------------------------
      ! Arrays with ghost cells
      !-------------------------------------------------
!     call alloc_x(q1,decomp_halo,.TRUE.)
!     call alloc_x(q2,decomp_halo,.TRUE.)
!     call alloc_x(q3,decomp_halo,.TRUE.)
!     call alloc_x(pr,decomp_halo,.TRUE.)
!     call alloc_x(dens,decomp_halo,.TRUE.)
!     call alloc_x(dph,decomp_halo,.TRUE.)
      allocate(q1(1:n3,xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1))
      allocate(q2(1:n3,xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1))
      allocate(q3(1:n3,xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1))
      allocate(pr(1:n3,xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1))
      allocate(dens(1:n3,xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1))
      allocate(dph(1:n3,xstart(2)-1:xend(2)+1,xstart(3)-1:xend(3)+1))

      !-----------------------------------------------
      ! Arrays without ghost cells
      !-----------------------------------------------
!     call alloc_x(rhs,ph,.TRUE.)
!     call alloc_x(ru1,ph,.TRUE.)
!     call alloc_x(ru2,ph,.TRUE.)
!     call alloc_x(ru3,ph,.TRUE.)
!     call alloc_x(ruro,ph,.TRUE.)
      allocate(rhs(1:n3,xstart(2):xend(2),xstart(3):xend(3)))
      allocate(dq(1:n3,xstart(2):xend(2),xstart(3):xend(3)))
      allocate(qcap(1:n3,xstart(2):xend(2),xstart(3):xend(3)))
      allocate(ru1(1:n3,xstart(2):xend(2),xstart(3):xend(3)))
      allocate(ru2(1:n3,xstart(2):xend(2),xstart(3):xend(3)))
      allocate(ru3(1:n3,xstart(2):xend(2),xstart(3):xend(3)))
      allocate(hro(1:n3,xstart(2):xend(2),xstart(3):xend(3)))
      allocate(ruro(1:n3,xstart(2):xend(2),xstart(3):xend(3)))
     
#ifdef STATS
      allocate(q1_me(1:n3m))
      allocate(q2_me(1:n3m))
      allocate(q3_me(1:n3m))
      allocate(q1_rms(1:n3m))
      allocate(q2_rms(1:n3m))
      allocate(q3_rms(1:n3m))
      allocate(dens_me(1:n3m))
      allocate(dens_rms(1:n3m))
      allocate(densq3_me(1:n3m))
#ifdef BALANCE
      allocate(disste(1:n3m))
      allocate(dissth(1:n3m))
#endif

#endif
          
      end subroutine mem_alloc

c==================================================      
      
      subroutine mem_dealloc
      use local_arrays
      use mpi_param
      use stat_arrays
      implicit none
      
      if(allocated(q1)) deallocate(q1)
      if(allocated(q2)) deallocate(q2)
      if(allocated(q3)) deallocate(q3)
      
      if(allocated(dens)) deallocate(dens)
      if(allocated(pr)) deallocate(pr)
      if(allocated(rhs)) deallocate(rhs)
      
      if(allocated(dph)) deallocate(dph)
      
      if(allocated(ru1)) deallocate(ru1)
      if(allocated(ru2)) deallocate(ru2)
      if(allocated(ru3)) deallocate(ru3)
      
      if(allocated(ruro)) deallocate(ruro)
      
      !---------------------------------------
      if(allocated(countj)) deallocate(countj)
      if(allocated(countk)) deallocate(countk)

      if(allocated(offsetj)) deallocate(offsetj)
      if(allocated(offsetk)) deallocate(offsetk)
      
      if(allocated(countf)) deallocate(countf)
      
      if(allocated(offsetf)) deallocate(offsetf)
      
#ifdef STATS
      if(allocated(q1_me)) deallocate(q1_me)
      if(allocated(q2_me)) deallocate(q2_me)
      if(allocated(q3_me)) deallocate(q3_me)
      if(allocated(q1_rms)) deallocate(q1_rms)
      if(allocated(q2_rms)) deallocate(q2_rms)
      if(allocated(q3_rms)) deallocate(q3_rms)
      if(allocated(dens_me)) deallocate(dens_me)
      if(allocated(dens_rms)) deallocate(dens_rms)
#endif
    
      end subroutine mem_dealloc
c================================================
#ifdef OLD
      
c================================================
      
#endif

