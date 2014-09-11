!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: hdnlro.F90                                     !
!    CONTAINS: subroutine hdnlro                          !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the temperature.                                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ExplicitTermsTemp
      use param
      use local_arrays, only: q2,q3,dens,q1,hro
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,kc,ic
      integer :: km,kp,jm,jp,im,ip
      real    :: h32,h33,udx2,udx1,h31
      real    :: udx1q,udx2q
      real    :: dq31,dq32

      udx1=dx1*0.25
      udx2=dx2*0.25
      udx1q=dx1q/pec
      udx2q=dx2q/pec

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,q1,q2,q3,n3m) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udx1) &
!$OMP   SHARED(udx2,udx1q,udx2q,udx3c,dens,hro) &
!$OMP   PRIVATE(ic,jc,kc,im,ip,km,kp,jm,jp) &
!$OMP   PRIVATE(h31,h32,h33,dq31,dq32)
      do ic=xstart(3),xend(3)
      im=ic-1
      ip=ic+1
      do jc=xstart(2),xend(2)
      jm=jc-1
      jp=jc+1
      do kc=2,n3m
      km=kc-1
      kp=kc+1
!
!
!    rho q1 term
!
!
!                d  rho q_t 
!             -----------
!                d   t      
!
      h31=((q1(km,jc,ip)+q1(kc,jc,ip))*(dens(kc,jc,ip)+dens(kc,jc,ic))- &
           (q1(km,jc,ic)+q1(kc,jc,ic))*(dens(kc,jc,ic)+dens(kc,jc,im)) &
          )*udx1
!
!
!    rho q2 term
!
!
!                d  rho q_r 
!             -----------
!                d   r      
!
      h32=((q2(kc,jp,ic)+q2(km,jp,ic))*(dens(kc,jp,ic)+dens(kc,jc,ic))- &
           (q2(kc,jc,ic)+q2(km,jc,ic))*(dens(kc,jc,ic)+dens(kc,jm,ic)) &
          )*udx2
!
!    rho q3 term
!
!
!                 d  rho q_x 
!                -----------
!                 d   x      
!
      h33=((q3(kp,jc,ic)+q3(kc,jc,ic))*(dens(kp,jc,ic)+dens(kc,jc,ic))- &
           (q3(kc,jc,ic)+q3(km,jc,ic))*(dens(kc,jc,ic)+dens(km,jc,ic)) &
          )*udx3c(kc)*0.25d0
!
!
!   11 second derivatives of dens
!
            dq31=(dens(kc,jc,ip) &
             -2.0*dens(kc,jc,ic) &
                 +dens(kc,jc,im))*udx1q
      
!
!   22 second derivatives of dens
!
            dq32=(dens(kc,jp,ic) &
             -2.0*dens(kc,jc,ic) &
                 +dens(kc,jm,ic))*udx2q
!
            hro(kc,jc,ic) = -(h31+h32+h33)+dq31+dq32
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
