      subroutine mkfftplans
      use param
      use fftw_params
      implicit none

!RO   Make fftw plans only once
!      call dfftw_plan_dft_r2c_1d(fwd_plan_y,m2m,yr,ya,FFTW_PATIENT)
!      call dfftw_plan_dft_c2r_1d(bck_plan_y,m2m,ya,yr,FFTW_PATIENT)
!
!      call dfftw_plan_dft_1d(fwd_plan_z,m1m,zr,za,FFTW_FORWARD,
!     %      FFTW_PATIENT)
!      call dfftw_plan_dft_1d(bck_plan_z,m1m,za,zr,FFTW_BACK,
!     %      FFTW_PATIENT)

      return
      end

