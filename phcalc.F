************************************************************************
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1to use the real fourier transform
c
      subroutine phcalc
      use param
      use fftw_params
      use local_arrays, only: dph
      use decomp_2d
      use decomp_2d_fft
      use mpih
      implicit none
      integer :: i,j,k,info
      complex :: acphT_b
      complex :: appph(m3m-2)
      complex :: amphT(m3m-1), apphT(m3m-1)
      complex, dimension(m3m) :: acphT,drhs,apph,amph
      integer :: phpiv(m3m)
      integer :: n2mh
      real,allocatable,dimension(:,:,:) :: ry1,rz1
      complex,allocatable,dimension(:,:,:) :: cy1,cz1,dphc
      logical,save:: printed=.false.

      integer ic,jc,kc

      integer,dimension(2)::dims,dummy_periods,dummy_coords
      integer ierror


!RO   Stuff for unoptimized FFTW


      allocate(ry1(ph%yst(1):ph%yen(1),
     &             ph%yst(2):ph%yen(2),
     &             ph%yst(3):ph%yen(3)))
      allocate(rz1(ph%zst(1):ph%zen(1),
     &             ph%zst(2):ph%zen(2),
     &             ph%zst(3):ph%zen(3)))
      allocate(cy1(sp%yst(1):sp%yen(1),
     &             sp%yst(2):sp%yen(2),
     &             sp%yst(3):sp%yen(3)))
      allocate(cz1(sp%zst(1):sp%zen(1),
     &             sp%zst(2):sp%zen(2),
     &             sp%zst(3):sp%zen(3)))
      allocate(dphc(sp%xst(1):sp%xen(1),
     &             sp%xst(2):sp%xen(2),
     &             sp%xst(3):sp%xen(3)))

      n2mh=n2m/2+1


      call transpose_x_to_y(dph(ph%xst(1):ph%xen(1),
     & ph%xst(2):ph%xen(2),
     & ph%xst(3):ph%xen(3)),ry1,ph)


!RO   Unoptimized FFT, needs to be written properly sometime

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(ph,n2mh)
!$OMP$ SHARED(fwd_plan_y,ry1,cy1)
!$OMP$ PRIVATE(ic,kc,yr,ya)
      do kc=ph%yst(1),ph%yen(1)
      do ic=ph%yst(3),ph%yen(3)
        yr(:)=ry1(kc,:,ic)
        call dfftw_execute_dft_r2c(fwd_plan_y,yr,ya)
        cy1(kc,1:n2mh,ic) = ya(1:n2mh)
       end do
      end do
!$OMP  END PARALLEL DO

      call transpose_y_to_z(cy1,cz1,sp)

!RO   Again shitty FFT, unoptimized


      do kc=sp%zst(1),sp%zen(1)
       do jc=sp%zst(2),sp%zen(2)
        zr(1:n1m)=cz1(kc,jc,1:n1m)
        call dfftw_execute(fwd_plan_z,zr,za)
        cz1(kc,jc,1:n1m) = za(1:n1m)
       end do
      end do



!EP   Normalize. FFT does not do this
      cz1 = cz1 / (n1m*n2m)

      if (.not.printed)then
        call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, 
     &    dims, dummy_periods, dummy_coords, ierror)
         printed=.true.
         do kc=0,dims(1)-1
         do ic=0,dims(2)-1
!JD no checks on x- or z-dimension, since we transform from z- into x-pencils, so these always overlap.
           if (sp%zst(2).le.sp%y1en(kc) .and. 
     &          sp%zen(2).ge.sp%y1st(kc)) then
!JD send to process with x-pencil defined by (kc,ic)
!JD x-bounds are taken from the z-pencils
!             send: sp%zst(1):sp%zen(1)
!JD y-bounds are the overlapping region of both pencils.
!                   max(sp%zst(2),sp%y1st(kc)):min(sp%zen(2),sp%y1en(kc))
!JD z-bounds are taken from the x-pencils.
!                   sp%z2st(ic):sp%z2en(ic)
             write(*,'(a,i4,i4)')'kc,ic=',kc,ic
             write(*,'(a,i4,i4)')'sp%zst(1):sp%zen(1)=',sp%zst(1),sp%zen(1)
             write(*,'(a,i4,i4)')'max(sp%zst(2),sp%y1st(kc)):min(sp%zen(2),sp%y1en(kc))=',
     &         max(sp%zst(2),sp%y1st(kc)),min(sp%zen(2),sp%y1en(kc))
             write(*,'(a,i4,i4)')'sp%z2st(ic):sp%z2en(ic)=',sp%z2st(ic),sp%z2en(ic)
             write(*,'(a,i4,i4,i4,i4,i4,i4)')'cz1(:,:,:)=',sp%zst(1),sp%zen(1),sp%zst(2),sp%zen(2),sp%zst(3),sp%zen(3)
           endif
!JD      Recvs:
!JD no checks on x- or z-dimension, since we transform from z- into x-pencils, so these always overlap.
           if (sp%xst(2).le.sp%y2en(ic) .and. 
     &          sp%xen(2).ge.sp%y2st(ic)) then
!JD recv from process with z-pencil defined by (kc,ic)
!JD x-bounds are taken from the z-pencils
!             send: sp%x1st(kc):sp%x1en(kc)
!JD y-bounds are the overlapping region of both pencils.
!                   max(sp%xst(2),sp%y2st(ic)):min(sp%xen(2),sp%y2en(ic))
!JD z-bounds are taken from the x-pencils.
!                   sp%xst(3):sp%xen(3)
             write(*,'(a,i4,i4)')'kc,ic=',kc,ic
             write(*,'(a,i4,i4)')'sp%x1st(kc):sp%x1en(kc)=',sp%x1st(kc),sp%x1en(kc)
             write(*,'(a,i4,i4)')'max(sp%xst(2),sp%y2st(ic)):min(sp%xen(2),sp%y2en(ic))=',
     &         max(sp%xst(2),sp%y2st(ic)),min(sp%xen(2),sp%y2en(ic))
             write(*,'(a,i4,i4)')'sp%xst(3):sp%xen(3)=',sp%xst(3),sp%xen(3)
             write(*,'(a,i4,i4,i4,i4,i4,i4)')'dphc(:,:,:)=',sp%xst(1),sp%xen(1),sp%xst(2),sp%xen(2),sp%xst(3),sp%xen(3)
           endif
         enddo
         enddo
      endif
      call transpose_z_to_y(cz1,cy1,sp)
      call transpose_y_to_x(cy1,dphc,sp)

!RO   Solve the tridiagonal matrix with complex coefficients

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(sp,n3m)
!$OMP$ SHARED(acphk,ak2,ak1,dphc,apphk,amphk)
!$OMP$ PRIVATE(drhs,apph,amph,acphT,acphT_b)
!$OMP$ PRIVATE(amphT,apphT,phpiv,info,appph)
      do i=sp%xst(3),sp%xen(3)
        do j=sp%xst(2),sp%xen(2)
         do k = 1,n3m
          acphT_b=1.0/(acphk(k)-ak2(j)-ak1(i))
          drhs(k)=dphc(k,j,i)*acphT_b
          apph(k)=apphk(k)*acphT_b
          amph(k)=amphk(k)*acphT_b
          acphT(k)=1.0d0
         enddo
  
         amphT=amph(2:n3m)
         apphT=apph(1:(n3m-1))

         call zgttrf(n3m, amphT, acphT, apphT, appph, phpiv, info)

         call zgttrs('N',n3m,1,amphT,acphT,apphT,appph,phpiv,drhs,
     &                 n3m, info)

          do k=1,n3m
            dphc(k,j,i) = drhs(k)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO


      call transpose_x_to_y(dphc,cy1,sp)

      call transpose_y_to_z(cy1,cz1,sp)

!RO   Inverse shitty FFT 1

      do kc=sp%zst(1),sp%zen(1)
      do jc=sp%zst(2),sp%zen(2)
        za(1:m1m)=cz1(kc,jc,1:m1m)
        call dfftw_execute(bck_plan_z,za,zr)
        cz1(kc,jc,1:m1m) = zr(1:m1m)
       end do
      end do

      call transpose_z_to_y(cz1,cy1,sp)


!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(bck_plan_y,ry1,cy1,sp)
!$OMP$ PRIVATE(ic,kc,yr,ya)
      do kc=sp%yst(1),sp%yen(1)
      do ic=sp%yst(3),sp%yen(3)
        ya(1:m2mh)=cy1(kc,sp%yst(2):sp%yen(2),ic)
        call dfftw_execute_dft_c2r(bck_plan_y,ya,yr)
        ry1(kc,1:m2m,ic) = yr(1:m2m)
       end do
      end do
!$OMP  END PARALLEL DO

      call transpose_y_to_x(ry1,dph(:,ph%xst(2):ph%xen(2),
     & ph%xst(3):ph%xen(3)),ph)


      if(allocated(dphc)) deallocate(dphc)
      if(allocated(rz1)) deallocate(rz1)
      if(allocated(cz1)) deallocate(cz1)
      if(allocated(ry1)) deallocate(ry1)
      if(allocated(cy1)) deallocate(cy1)


      return
      end
