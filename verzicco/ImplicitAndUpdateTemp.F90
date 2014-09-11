!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: invtrro.F90                                    !
!    CONTAINS: subroutine invtrro                         !
!                                                         ! 
!    PURPOSE: Compute the linear terms associated to      !
!     the temperature and calls the implicit solver.      !                                 !
!     After this routine, the temperature field has been  ! 
!     updated to the new timestep                         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ImplicitAndUpdateTemp
      use param
      use local_arrays, only: dens,hro,ruro,rhs
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,kc,ic
      integer :: km,kp
      real    :: alpec,dq33
      real    :: app,acc,amm



      alpec=al/pec

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,n3m,dens) &
!$OMP   SHARED(kmv,kpv,am3ck,ac3ck,ap3ck) &
!$OMP   SHARED(ga,ro,alpec,dt) &
!$OMP   SHARED(rhs,ruro,hro) &
!$OMP   PRIVATE(ic,jc,kc,km,kp) &
!$OMP   PRIVATE(amm,acc,app) &
!$OMP   PRIVATE(dq33)
      do ic=xstart(3),xend(3)
      do jc=xstart(2),xend(2)
      do kc=2,n3m

!   Calculate second derivative of temperature in the z-direction.
!   This is the only term calculated implicitly for temperature.

               dq33= dens(kc+1,jc,ic)*ap3ck(kc) &
                    +dens(kc  ,jc,ic)*ac3ck(kc) &
                    +dens(kc-1,jc,ic)*am3ck(kc)


!    Calculate right hand side of Eq. 5 (VO96)

            rhs(kc,jc,ic)=(ga*hro(kc,jc,ic)+ro*ruro(kc,jc,ic) &
                    +alpec*dq33)*dt

!    Store the non-linear terms for the calculation of 
!    the next timestep

            ruro(kc,jc,ic)=hro(kc,jc,ic)

        enddo
       enddo
      enddo
!$OMP END PARALLEL DO


!  Solve equation and update temperature

      call SolveImpEqnUpdate_Temp

!  Set boundary conditions on the temperature field at top
!  and bottom plates. This seems necessary.

       dens(1,xstart(2):xend(2),xstart(3):xend(3)) &
          = denbs(xstart(2):xend(2),xstart(3):xend(3))

       dens(n3,xstart(2):xend(2),xstart(3):xend(3)) &
          = denbn(xstart(2):xend(2),xstart(3):xend(3))


      return
      end
