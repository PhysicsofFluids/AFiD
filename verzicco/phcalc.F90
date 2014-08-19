!************************************************************************
!  this subroutine perform the calculation of dph , periodic direction
!  along x3 and x1to use the real fourier transform
!
      subroutine phcalc
      use, intrinsic :: iso_c_binding
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

      type(fftw_iodim),dimension(1) :: iodim
      type(fftw_iodim),dimension(2) :: iodim_howmany

!RO   Stuff for unoptimized FFTW


      allocate(ry1(ph%yst(1):ph%yen(1),                                 &
     &             ph%yst(2):ph%yen(2),                                 &
     &             ph%yst(3):ph%yen(3)))
      allocate(rz1(ph%zst(1):ph%zen(1),                                 &
     &             ph%zst(2):ph%zen(2),                                 &
     &             ph%zst(3):ph%zen(3)))
      allocate(cy1(sp%yst(1):sp%yen(1),                                 &
     &             sp%yst(2):sp%yen(2),                                 &
     &             sp%yst(3):sp%yen(3)))
      allocate(cz1(sp%zst(1):sp%zen(1),                                 &
     &             sp%zst(2):sp%zen(2),                                 &
     &             sp%zst(3):sp%zen(3)))
      allocate(dphc(sp%xst(1):sp%xen(1),                                &
     &             sp%xst(2):sp%xen(2),                                 &
     &             sp%xst(3):sp%xen(3)))

      n2mh=n2m/2+1

!EP   Here array temporaries are created, this should be avoided
!     call transpose_x_to_y(dph(ph%xst(1):ph%xen(1),
!    & ph%xst(2):ph%xen(2),
!    & ph%xst(3):ph%xen(3)),ry1,ph)

      call transpose_x_to_y(dph,ry1,ph)

      if (.not.planned) then
        iodim(1)%n=n1m
        iodim(1)%is=(sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)
        iodim(1)%os=(sp%zen(1)-sp%zst(1)+1)*(sp%zen(2)-sp%zst(2)+1)
        iodim_howmany(1)%n=(sp%zen(1)-sp%zst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(sp%zen(2)-sp%zst(2)+1)
        iodim_howmany(2)%is=(sp%zen(1)-sp%zst(1)+1)
        iodim_howmany(2)%os=(sp%zen(1)-sp%zst(1)+1)
        fwd_guruplan_z=fftw_plan_guru_dft(1,iodim,                      &
     &    2,iodim_howmany,cz1,cz1,                                      &
     &    FFTW_FORWARD,FFTW_PATIENT)
        iodim(1)%n=m1m
        bwd_guruplan_z=fftw_plan_guru_dft(1,iodim,                      &
     &    2,iodim_howmany,cz1,cz1,                                      &
     &    FFTW_BACKWARD,FFTW_PATIENT)

        iodim(1)%n=m2m
        iodim(1)%is=ph%yen(1)-ph%yst(1)+1
        iodim(1)%os=sp%yen(1)-sp%yst(1)+1
        iodim_howmany(1)%n=(ph%yen(1)-ph%yst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(ph%yen(3)-ph%yst(3)+1)
        iodim_howmany(2)%is=(ph%yen(1)-ph%yst(1)+1)                     &
     &    *(ph%yen(2)-ph%yst(2)+1)
        iodim_howmany(2)%os=(sp%yen(1)-sp%yst(1)+1)                     &
     &    *(sp%yen(2)-sp%yst(2)+1)
        fwd_guruplan_y=fftw_plan_guru_dft_r2c(1,iodim,                  &
     &    2,iodim_howmany,ry1,cy1,                                      &
     &    FFTW_PATIENT)

        iodim(1)%n=m2m
        iodim(1)%is=sp%yen(1)-sp%yst(1)+1
        iodim(1)%os=ph%yen(1)-ph%yst(1)+1
        iodim_howmany(1)%n=(sp%yen(1)-sp%yst(1)+1)
        iodim_howmany(1)%is=1
        iodim_howmany(1)%os=1
        iodim_howmany(2)%n=(sp%yen(3)-sp%yst(3)+1)
        iodim_howmany(2)%is=(sp%yen(1)-sp%yst(1)+1)                     &
     &    *(sp%yen(2)-sp%yst(2)+1)
        iodim_howmany(2)%os=(ph%yen(1)-ph%yst(1)+1)                     &
     &    *(ph%yen(2)-ph%yst(2)+1)
        bwd_guruplan_y=fftw_plan_guru_dft_c2r(1,iodim,                  &
     &    2,iodim_howmany,cy1,ry1,                                      &
     &    FFTW_PATIENT)
        planned=.true.
      endif

      call dfftw_execute_dft_r2c(fwd_guruplan_y,ry1,cy1)

      call transpose_y_to_z(cy1,cz1,sp)

      call dfftw_execute_dft(fwd_guruplan_z,cz1,cz1)


!EP   Normalize. FFT does not do this
      cz1 = cz1 / (n1m*n2m)

      call transpose_z_to_x(cz1,dphc,sp)

!RO   Solve the tridiagonal matrix with complex coefficients

!$OMP  PARALLEL DO                                                      &
!$OMP   DEFAULT(none)                                                   &
!$OMP   SHARED(sp,n3m)                                                  &
!$OMP   SHARED(acphk,ak2,ak1,dphc,apphk,amphk)                          &
!$OMP   PRIVATE(drhs,apph,amph,acphT,acphT_b)                           &
!$OMP   PRIVATE(amphT,apphT,phpiv,info,appph)
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

         call zgttrs('N',n3m,1,amphT,acphT,apphT,appph,phpiv,drhs,      &
     &                 n3m, info)

          do k=1,n3m
            dphc(k,j,i) = drhs(k)
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO

      call transpose_x_to_z(dphc,cz1,sp)

      call dfftw_execute_dft(bwd_guruplan_z,cz1,cz1)

      call transpose_z_to_y(cz1,cy1,sp)

      call dfftw_execute_dft_c2r(bwd_guruplan_y,cy1,ry1)

!EP   Here array temporaries are created, this should be avoided
!     call transpose_y_to_x(ry1,dph(:,ph%xst(2):ph%xen(2),
!    & ph%xst(3):ph%xen(3)),ph)
      call transpose_y_to_x(ry1,dph,ph)


      if(allocated(dphc)) deallocate(dphc)
      if(allocated(rz1)) deallocate(rz1)
      if(allocated(cz1)) deallocate(cz1)
      if(allocated(ry1)) deallocate(ry1)
      if(allocated(cy1)) deallocate(cy1)


      return
      end
