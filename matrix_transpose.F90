!===========================================
      subroutine PackZ_UnpackR(aa,bb)
      use mpih
      use mpi_param
      use param, only: n2m,n3m,n2,n3,n1
      implicit none
      real, intent(in) :: aa(1:n1,1:n2,kstart:kend)
      real, intent(out) :: bb(1:n3,1:n1,jstart:jend)
      real, allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:)
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countj(i)*n1
        aak(i)= dj* countk(i)*n1
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
     
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1*n2m*dk-1),stat=merr)
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(rbuf)) allocate(rbuf(0:n1*n3m*dj-1),stat=merr)
      
      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

      nc=0
      do kk = 0, numtasks-1
        dr= countj(kk)
        offsetr = offsetj(kk)
        do k = kstart,kend
          do j=1,dr
            do i=1,n1
              sbuf(nc) = aa(i,j+offsetr,k)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo
      call MPI_ALLTOALLV(sbuf, aaj,dispj, MDP,
     &                  rbuf, aak,dispk, MDP,
     &                  MPI_COMM_WORLD, ierr)
     
     
      nc=0
      do kk = 0, numtasks-1
        dz= countk(kk)
        offsetz = offsetk(kk)
        do k = 1,dz
          do j=jstart,jend
            do i=1,n1
              bb(k+offsetz,i,j) = rbuf(nc)
              nc=nc+1
            enddo
          enddo
        enddo
      enddo

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackZ_UnpackR
 
!============================================
      subroutine PackR_UnpackZ(aa,bb)
      use mpih
      use mpi_param
      use param, only: n2m,n3m,n2,n3,n1
      implicit none
      real, intent(in) :: aa(1:n3,1:n1,jstart:jend)
      real, intent(out) :: bb(1:n1,1:n2,kstart:kend)
      real, allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countj(i)*n1
        aak(i)= dj* countk(i)*n1
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(.not. allocated(rbuf)) allocate(rbuf(0:n1*n2m*dk-1),stat=merr)

      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1*n3m*dj-1),stat=merr)

      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory for rbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 

      nc=0
      do kk = 0, numtasks-1
        dz= countk(kk)
        offsetz = offsetk(kk)
          do k = 1,dz
            do j=jstart,jend
              do i=1,n1
              sbuf(nc) = aa(k+offsetz,i,j) 
              nc=nc+1
              enddo
          enddo
        enddo
      enddo
      
      call MPI_ALLTOALLV(sbuf, aak,dispk, MDP,
     &                  rbuf, aaj,dispj, MDP,
     &                  MPI_COMM_WORLD, ierr)
     
      nc=0
      do kk = 0, numtasks-1
        dr= countj(kk)
        offsetr = offsetj(kk)
          do k = kstart,kend
            do j=1,dr
              do i=1,n1
              bb(i,j+offsetr,k) = rbuf(nc)
              nc=nc+1
              enddo
          enddo
        enddo
      enddo

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackR_UnpackZ
!==========================================
      subroutine PackZ_UnpackRP(aa,bb)
      use mpih
      use mpi_param
      use param, only: n2,n3m,n2,n3,n1
      implicit none
      real,intent(in) :: aa(1:n1,1:n2+1,kstart:kend)
      real,intent(out) :: bb(1:n3,1:n1,jstartp:jendp)
      real, allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countjp(i)*n1
        aak(i)= djp* countk(i)*n1
      end do
      
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 

      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
 
     
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1*(n2+1)*dk-1)
     & ,stat=merr)

      if(.not. allocated(rbuf)) allocate(rbuf(0:n1*n3m*djp-1)
     & ,stat=merr)

      nc=0
      do kk = 0, numtasks-1
        dr= countjp(kk)
        offsetr = offsetjp(kk)
        do k = kstart,kend
          do j=1,dr
              do i=1,n1
              sbuf(nc) = aa(i,j+offsetr,k)
              nc=nc+1
              enddo
          enddo
        enddo
      enddo
      
      call MPI_ALLTOALLV(sbuf, aaj,dispj, MDP,
     &                  rbuf, aak,dispk, MDP,
     &                  MPI_COMM_WORLD, ierr)
     
     
      nc=0
      do kk = 0, numtasks-1
        dz= countk(kk)
        offsetz = offsetk(kk)
        do k = 1,dz
          do j=jstartp,jendp
              do i=1,n1
              bb(k+offsetz,i,j) = rbuf(nc)
              nc=nc+1
              enddo
          enddo
        enddo
      enddo

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackZ_UnpackRP
 
!============================================
      subroutine PackR_UnpackZP(aa,bb)
      use mpih
      use mpi_param
      use param, only: n2,n3m,n3,n1
      implicit none
      real, intent(in) :: aa(1:n3,1:n1,jstartp:jendp)
      real, intent(out) :: bb(1:n1,1:n2+1,kstart:kend)
      real, allocatable :: sbuf(:),rbuf(:)
      integer, allocatable :: aaj(:), aak(:)
      integer, allocatable :: dispj(:), dispk(:) 
      integer :: dr,dz,offsetr,offsetz
      integer :: i,j,k,kk,nc
      integer :: merr

      allocate(aaj(0:numtasks-1))
      allocate(aak(0:numtasks-1))

      do i=0,numtasks-1
        aaj(i)= dk* countjp(i)*n1
        aak(i)= djp* countk(i)*n1
      end do
       
      allocate(dispj(0:numtasks-1))
      allocate(dispk(0:numtasks-1)) 
      
      dispj(:)=0
      dispk(:)=0
      do i=1,numtasks-1
        dispj(i)= dispj(i-1) + aaj(i-1)
        dispk(i)= dispk(i-1) + aak(i-1)
      end do
      
      if(.not. allocated(rbuf)) allocate(rbuf(0:n1*(n2+1)*dk-1)
     & ,stat=merr)

      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory for sbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      if(.not. allocated(sbuf)) allocate(sbuf(0:n1*n3m*djp-1),
     & stat=merr)

      if(merr .ne. 0) then
        write(6,*)"process  ",myid," failed to allocate memory for rbuf"
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif 
      
      nc=0
      do kk = 0, numtasks-1
        dz= countk(kk)
        offsetz = offsetk(kk)
          do k = 1,dz
            do j=jstartp,jendp
              do i=1,n1
              sbuf(nc) = aa(k+offsetz,i,j) 
              nc=nc+1
              enddo
          enddo
        enddo
      enddo
      
      
      call MPI_ALLTOALLV(sbuf, aak,dispk, MDP,
     &                  rbuf, aaj,dispj, MDP,
     &                  MPI_COMM_WORLD, ierr)
     
      nc=0
      do kk = 0, numtasks-1
        dr= countjp(kk)
        offsetr = offsetjp(kk)
          do k = kstart,kend
            do j=1,dr
              do i=1,n1
              bb(i,j+offsetr,k) = rbuf(nc)
              nc=nc+1
              enddo
          enddo
        enddo
      enddo

      if(allocated(sbuf)) deallocate(sbuf)
      if(allocated(rbuf)) deallocate(rbuf)
     
      if(allocated(aaj)) deallocate(aaj)
      if(allocated(aak)) deallocate(aak)
     
      if(allocated(dispj)) deallocate(dispj)
      if(allocated(dispk)) deallocate(dispk)
      
      end subroutine PackR_UnpackZP
!==========================================
