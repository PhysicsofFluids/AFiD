!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolvePressureCorrection.F90                    !
!    CONTAINS: subroutine SolvePressureCorrection,        !
!     CreateFFTTmpArrays, DestroyFFTTmpArrays             ! 
!                                                         ! 
!    PURPOSE: Compute the pressure correction by solving  !
!     a Poisson equation                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SolvePressureCorrection
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
      complex :: appph(nxm-2)
      complex, dimension(nxm) :: acphT,apph,amph
      integer :: phpiv(nxm)
      integer :: nymh

      type(fftw_iodim),dimension(1) :: iodim
      type(fftw_iodim),dimension(2) :: iodim_howmany

!RO   Allocate variables for FFT transform

      call CreateFFTTmpArrays

      nymh=nym/2+1

      call transpose_x_to_y(dph,ry1,ph)

!RO   Plan FFT transforms if not planned previously 

      if (.not.planned) then
        iodim(1)%n=nzm
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
     &    FFTW_FORWARD,FFTW_ESTIMATE)
        iodim(1)%n=nzm
        bwd_guruplan_z=fftw_plan_guru_dft(1,iodim,                      &
     &    2,iodim_howmany,cz1,cz1,                                      &
     &    FFTW_BACKWARD,FFTW_ESTIMATE)

        if (.not.c_associated(bwd_guruplan_z)) then
          if (ismaster) print*,'Failed to create guru plan. You should'
          if (ismaster) print*,'link with FFTW3 before MKL'
          if (ismaster) print*,'Please check linking order.'
          call MPI_Abort(MPI_COMM_WORLD,1,info)
        endif

        iodim(1)%n=nym
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
     &    FFTW_ESTIMATE)

        iodim(1)%n=nym
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
     &    FFTW_ESTIMATE)
        planned=.true.
      endif

      call dfftw_execute_dft_r2c(fwd_guruplan_y,ry1,cy1)

      call transpose_y_to_z(cy1,cz1,sp)

      call dfftw_execute_dft(fwd_guruplan_z,cz1,cz1)

!EP   Normalize. FFT does not do this
      cz1 = cz1 / (nzm*nym)

      call transpose_z_to_x(cz1,dphc,sp)

!RO   Solve the tridiagonal matrix with complex coefficients

!$OMP  PARALLEL DO COLLAPSE(2)                                          &
!$OMP   DEFAULT(none)                                                   &
!$OMP   SHARED(sp,nxm)                                                  &
!$OMP   SHARED(acphk,ak2,ak1,dphc,apphk,amphk)                          &
!$OMP   PRIVATE(apph,amph,acphT,acphT_b)                                &
!$OMP   PRIVATE(phpiv,info,appph)
      do i=sp%xst(3),sp%xen(3)
        do j=sp%xst(2),sp%xen(2)
         do k = 1,nxm
          acphT_b=1.0/(acphk(k)-ak2(j)-ak1(i))
          dphc(k,j,i)=dphc(k,j,i)*acphT_b
          apph(k)=apphk(k)*acphT_b
          amph(k)=amphk(k)*acphT_b
          acphT(k)=1.0d0
         enddo
  

         call zgttrf(nxm, amph(2), acphT, apph(1), appph, phpiv, info)

         if (info.gt.0) then
           print*,'Singular value found in LAPACK routine zgttrf: info=',info
           print*,'Please try to adjust either NX or STR3 in bou.in'
           call MPI_Abort(MPI_COMM_WORLD,1,ierr)
         endif

         call zgttrs('N',nxm,1,amph(2),acphT,apph(1),appph,phpiv,      &
                       dphc(1,j,i), nxm, info)

        enddo
      enddo
!$OMP END PARALLEL DO

      call transpose_x_to_z(dphc,cz1,sp)

      call dfftw_execute_dft(bwd_guruplan_z,cz1,cz1)

      call transpose_z_to_y(cz1,cy1,sp)

      call dfftw_execute_dft_c2r(bwd_guruplan_y,cy1,ry1)

      call transpose_y_to_x(ry1,dph,ph)

      call DestroyFFTTmpArrays

      return
      end

!======================================================================

      subroutine CreateFFTTmpArrays
      use fftw_params
      use decomp_2d
      use decomp_2d_fft
      implicit none

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

      return
      end

!======================================================================

      subroutine DestroyFFTTmpArrays
      use fftw_params
      implicit none

      if(allocated(dphc)) deallocate(dphc)
      if(allocated(rz1)) deallocate(rz1)
      if(allocated(cz1)) deallocate(cz1)
      if(allocated(ry1)) deallocate(ry1)
      if(allocated(cy1)) deallocate(cy1)

      return
      end
