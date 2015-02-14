      subroutine SlabDumper
      use param
      use local_arrays, only: temp,vz,vy,vx
      use stat3_param
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: i,j,m
      real,dimension(xstart(2):xend(2),xstart(3):xend(3)) :: &
     &      vxcc,vzcc,vycc,tempcc
      character*70 :: filnam
      character*1 :: charm

!EP   Slabs
!EP   cell center only vx

      do m=1,9
!$OMP  PARALLEL DO DEFAULT(SHARED) &
!$OMP   PRIVATE(i,j)
        do i=xstart(3),xend(3)
         do j=xstart(2),xend(2)
           vxcc(j,i) = (vx(kslab(m),j,i)+vx(kslab(m)+1,j,i))*0.5
           vycc(j,i) = vy(kslab(m),j,i)
           vzcc(j,i) = vz(kslab(m),j,i)
           tempcc(j,i) = temp(kslab(m),j,i)
          enddo
         enddo
!$OMP  END PARALLEL DO
      write(charm,28) m
   28 format(i1.1)
      filnam='slab'//charm//'vx_'
      call DumpSingleSlab(vxcc,filnam)
      filnam='slab'//charm//'vy_'
      call DumpSingleSlab(vycc,filnam)
      filnam='slab'//charm//'vz_'
      call DumpSingleSlab(vzcc,filnam)
      filnam='slab'//charm//'temp_'
      call DumpSingleSlab(tempcc,filnam)
      enddo

      return
      end subroutine SlabDumper

!===========================================================================

      subroutine InitializeSlabDump
      use param
      use stat3_param
      implicit none
      integer :: i,k,j
      real :: xmloc
      character(len=4) :: dummy

!EP   Read from stst3.in
      
      open(unit=19,file='stst3.in',status='old')
        read(19,301) dummy
        read(19,*) (zslab(i),i=2,9)
301     format(a4)                
      close(19)

!EP   Compute which kslab corresponds to which zslab
      
      kslab = 2
      
        do k=2,nxm
          xmloc=xm(k)
          do j=2,9
            if(xm(k).gt.zslab(j).and.xm(k-1).lt.zslab(j)) then
             kslab(j) = k
            endif
          enddo
        enddo


!EP   Write probe and slab locations
      
      if (ismaster) then
      open(unit=23,file='stst3locs.out',status='unknown')
        rewind(23)
        write(23,*) (kslab(i),i=1,9)
      close(23)
      endif

      return
      end subroutine InitializeSlabDump

!==================================================================
      
      subroutine DumpSingleSlab(var,filnam)
      USE param
      use mpih
      USE hdf5
      use decomp_2d, only: xstart,xend
      IMPLICIT none

      real, intent(in) :: var(xstart(2):xend(2) &
     &                  ,xstart(3):xend(3))

      real :: tprfi
      integer :: itime

      character*70,intent(in) :: filnam
      character*70 :: namfile,dsetname
      character*8 :: ipfi

      tprfi = 1/tout
      itime=nint(time*tprfi)
      write(ipfi,82)itime
   82 format(i8.8)

      namfile=trim('./stst3/'//trim(filnam)//trim(ipfi)//'.h5')
      dsetname = trim('var')

      call HdfWriteReal2D(dsetname,namfile,var)

      if(ismaster) then
       dsetname = trim('time')
       call HdfSerialWriteRealScalar(dsetname,namfile,time)
      endif

      return                                                          
      end subroutine DumpSingleSlab
