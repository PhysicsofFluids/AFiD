!***********************************************************************
      subroutine inirea
      use mpih
      use decomp_2d
      use local_arrays
      use param
      IMPLICIT NONE
      real ::  aaa
      character*70 :: filcnw2
      integer :: ihist
      integer :: n1om,n2om,n3o,n2o,n1o
      integer :: xs2og,xe2og,xs3og,xe3og
      integer :: istro3
      integer (kind=MPI_ADDRESS_KIND) :: extent,lb
      real :: stro3
      real :: intinfo(1:4)
      real, allocatable, dimension(:,:,:) :: densold,q2old,q3old,q1old
      
!EP   Reading old grid information by rank0
      if (nrank .eq. 0) then
      filcnw2 = 'continua_grid.dat'
      open(13,file=filcnw2,status='unknown')
      rewind(13)                                                      
      read(13,*) n1o,n2o,n3o
      read(13,*) aaa,time
      read(13,*) istro3,stro3
      close(13)
      endif
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!EP   Bcasting old grid information and time
      call MPI_BCAST(n1o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(n2o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(n3o,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(istro3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stro3,1,MDP,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(time,1,MDP,0,MPI_COMM_WORLD,ierr)
      
!EP   Check whether grid specifications have been updated
      if(n2o.ne.n2.or.n3o.ne.n3.or.n1o.ne.n1 &
       .or.istro3.ne.istr3.or.stro3.ne.str3) then
      if(nrank.eq.0) write(*,*) "Interpolating new grid"
      if(n1.gt.n1o*2.or.n2.gt.n2o*2.or.n3.gt.n3o*2) then
      if(nrank.eq.0) write(*,*) "New grid resolution cannot be more ", &
       "than twice the old resolution"
      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      endif

      n1om = n1o - 1
      n2om = n2o - 1
      
      intinfo(1) = istro3
      intinfo(2) = stro3

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      xs2og = floor(real(xstart(2)*n2om/n2m))
      xe2og = ceiling(real(xend(2)*n2om/n2m))
      xs3og = floor(real(xstart(3)*n1om/n1m))
      xe3og = ceiling(real(xend(3)*n1om/n1m))

      xs2og = max(xs2og,1)
      xe2og   = min(xe2og,n2om)
      xs3og = max(xs3og,1)
      xe3og   = min(xe3og,n1om)
      
      lb=0
      call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,extent,ierr)
!EP   dens
      allocate(densold(0:n3o+1,xs2og-1:xe2og+1,xs3og-1:xe3og+1))
      
      call hdf_read(n1o,n2o,n3o,xs2og,xe2og, &
       xs3og,xe3og,4,densold(1:n3o,xs2og-1:xe2og+1,xs3og-1:xe3og+1))

      call interp(densold,dens(1:n3,xstart(2):xend(2),xstart(3):xend(3)) &
       ,n1o,n2o,n3o,istro3,stro3,4,xs2og,xe2og,xs3og,xe3og)

      deallocate(densold)

!EP   q1
      allocate(q1old(0:n3o+1,xs2og-1:xe2og+1,xs3og-1:xe3og+1))
      
      call hdf_read(n1o,n2o,n3o,xs2og,xe2og, &
       xs3og,xe3og,1,q1old(1:n3o,xs2og-1:xe2og+1,xs3og-1:xe3og+1))

      call interp(q1old,q1(1:n3,xstart(2):xend(2),xstart(3):xend(3)) &
       ,n1o,n2o,n3o,istro3,stro3,1,xs2og,xe2og,xs3og,xe3og)

      deallocate(q1old)

!EP   q2
      allocate(q2old(0:n3o+1,xs2og-1:xe2og+1,xs3og-1:xe3og+1))
      
      call hdf_read(n1o,n2o,n3o,xs2og,xe2og, &
     & xs3og,xe3og,2,q2old(1:n3o,xs2og-1:xe2og+1,xs3og-1:xe3og+1))

      call interp(q2old,q2(1:n3,xstart(2):xend(2),xstart(3):xend(3)) &
     & ,n1o,n2o,n3o,istro3,stro3,2,xs2og,xe2og,xs3og,xe3og)

      deallocate(q2old)

!EP   q3
      allocate(q3old(0:n3o+1,xs2og-1:xe2og+1,xs3og-1:xe3og+1))
      
      call hdf_read(n1o,n2o,n3o,xs2og,xe2og, &
     & xs3og,xe3og,3,q3old(1:n3o,xs2og-1:xe2og+1,xs3og-1:xe3og+1))

      call interp(q3old,q3(1:n3,xstart(2):xend(2),xstart(3):xend(3)) &
     & ,n1o,n2o,n3o,istro3,stro3,3,xs2og,xe2og,xs3og,xe3og)

      deallocate(q3old)

      else

!EP   One to one HDF read
      call hdf_read(n1,n2,n3,xstart(2),xend(2) &
     & ,xstart(3),xend(3),4,dens)
      call hdf_read(n1,n2,n3,xstart(2),xend(2) &
     & ,xstart(3),xend(3),1,q1)
      call hdf_read(n1,n2,n3,xstart(2),xend(2) &
     & ,xstart(3),xend(3),2,q2)
      call hdf_read(n1,n2,n3,xstart(2),xend(2) &
     & ,xstart(3),xend(3),3,q3)

      endif

      if (ireset.eq.1) then                                             
       ihist=0                                                          
      time=0.
      endif                                                             

      return                                                            
      end                                                               
!                                                                       
