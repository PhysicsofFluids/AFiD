!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: TimeMarcher.F90                                !
!    CONTAINS: subroutine TimeMarcher                     !
!                                                         ! 
!    PURPOSE: Main time integrating routine, which calls  !
!     other subroutines for calculating the Navier-Stokes !
!     equations and advancing velocity and temperature in !
!     time                                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine TimeMarcher
      use param
      use local_arrays
      use mpih
      use decomp_2d
      implicit none
      integer :: ns
      integer :: j,k,i

      beta=dt/ren*0.5d0

      do ns=1,nsst                                                 

!RO     Coefficients for time marching integration (alpha, gamma, rho)

        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        call ExplicitTermsVX
        call ExplicitTermsVY
        call ExplicitTermsVZ
        call ExplicitTermsTemp     

        call ImplicitAndUpdateVX
        call ImplicitAndUpdateVY
        call ImplicitAndUpdateVZ
        call ImplicitAndUpdateTemp

        call update_halo(vy,lvlhalo)
        call update_halo(vz,lvlhalo)


        call CalcLocalDivergence
        call SolvePressureCorrection

!EP this copy can be avoided by changing transpose_x_to_y_real and
!transpose_y_to_x_real so these routines can handle arrays with
!halo. This copy is a defacto array temporary. Using inferred size
!arrays in the transpose calls results in 5 more of these, and more
!memory usage.  Time spent on this copy is 0.1% for 65^3 grid.

        do i=xstart(3),xend(3)
          do j=xstart(2),xend(2)
            do k=1,nxm
              dphhalo(k,j,i) = dph(k,j,i)
            enddo
          enddo
        enddo

        call update_halo(dphhalo,lvlhalo)

        call CorrectVelocity
        call CorrectPressure

        call update_halo(vx,lvlhalo)
        call update_halo(vy,lvlhalo)
        call update_halo(vz,lvlhalo)
        call update_halo(pr,lvlhalo)
        call update_halo(temp,lvlhalo)

        enddo


      return                                                            
      end                                                               
