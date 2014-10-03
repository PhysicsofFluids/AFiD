!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: tsch.F90                                       !
!    CONTAINS: subroutine tschem                          !
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
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        call ExplicitTermsVX
        call ExplicitTermsVY
        call ExplicitTermsVZ
        call ExplicitTermsTemp     

!        call Checksum3DArray(q2,cksum)
!        if(ismaster) write(6,*) 'cksum',cksum
!        call Checksum3DArray(q3,cksum)
!        if(ismaster) write(6,*) 'cksum',cksum
!        call Checksum3DArray(dens,cksum)
!        if(ismaster) write(6,*) 'cksum',cksum

        call ImplicitAndUpdateVX
        call ImplicitAndUpdateVY
        call ImplicitAndUpdateVZ
        call ImplicitAndUpdateTemp

        call update_halo(q1,1)
        call update_halo(q2,1)

        call CalculateLocalDivergence
        call SolvePressureCorrection

!EP this copy can be avoided by changing transpose_x_to_y_real and
!transpose_y_to_x_real so these routines can handles arrays with
!halo. This copy is a defacto array temporary. Using inferred size
!arrays in the transpose calls results in 5 more of these. Time spent on
!this copy is 0.1% for 65^3 grid.
        do i=xstart(3),xend(3)
          do j=xstart(2),xend(2)
            do k=1,nxm
              dphhalo(k,j,i) = dph(k,j,i)
            enddo
          enddo
        enddo

        call update_halo(dphhalo,1)

        call CorrectVelocity
        call CorrectPressure

        call update_halo(q1,1)
        call update_halo(q2,1)
        call update_halo(q3,1)
        call update_halo(pr,1)
        call update_halo(dens,1)

        enddo


      return                                                            
      end                                                               
