!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: hdnl3.F90                                      !
!    CONTAINS: subroutine hdnl3                           !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the vertical dimension (x3)         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine ExplicitTermsVX
      use param
      use local_arrays, only: vz,vy,vx,temp,qcap
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: jc,kc
      integer :: km,kp,jmm,jpp,ic,imm,ipp
      real    :: h32,h33,h31
      real    :: udz,udy,tempit
      real    :: udzq,udyq
      real    :: dq31,dq32

      udy=dy*0.25
      udz=dz*0.25

      udyq=dyq/ren
      udzq=dzq/ren

!$OMP  PARALLEL DO &
!$OMP   DEFAULT(none) &
!$OMP   SHARED(xstart,xend,nxm,vz,vy,vx,dz,dy) &
!$OMP   SHARED(kmv,kpv,am3sk,ac3sk,ap3sk,udz) &
!$OMP   SHARED(udy,udzq,udyq,udx3c,qcap,temp) &
!$OMP   PRIVATE(ic,jc,kc,imm,ipp,km,kp) &
!$OMP   PRIVATE(jmm,jpp,tempit) &
!$OMP   PRIVATE(h31,h32,h33,dq31,dq32)
      do ic=xstart(3),xend(3)
       imm=ic-1
       ipp=ic+1
       do jc=xstart(2),xend(2)
        jmm=jc-1
        jpp=jc+1
        do kc=2,nxm
         km=kc-1
         kp=kc+1
!
!    vx vz term
!
!
!                d  q_x q_t 
!             -----------
!                d   t      
!
!
      h31=(((vz(kc,jc,ipp)+vz(km,jc,ipp)) &
           *(vx(kc,jc,ipp)+vx(kc,jc,ic))) &
          -((vz(kc,jc,ic)+vz(km,jc,ic)) &
           *(vx(kc,jc,ic)+vx(kc,jc,imm))))*udz
!
!    vx vy term
!
!
!                d  q_x q_r 
!             -----------
!                d   r      
!
      h32=(((vy(kc,jpp,ic)+vy(km,jpp,ic)) &
           *(vx(kc,jpp,ic)+vx(kc,jc,ic))) &
          -((vy(kc,jc,ic)+vy(km,jc,ic)) &
           *(vx(kc,jc,ic)+vx(kc,jmm,ic))))*udy
!
!    vx vx term
!
!
!                 d  q_x q_x 
!                -----------
!                 d   x      
!
      h33=((vx(kp,jc,ic)+vx(kc,jc,ic))*(vx(kp,jc,ic)+vx(kc,jc,ic)) &
          -(vx(kc,jc,ic)+vx(km,jc,ic))*(vx(kc,jc,ic)+vx(km,jc,ic)) &
          )*udx3c(kc)*0.25d0
!
!  add the buoyancy term
!
          tempit=temp(kc,jc,ic)

!
!   11 second derivatives of vx
!
            dq31=(vx(kc,jc,imm) &
                 -2.0*vx(kc,jc,ic) &
                 +vx(kc,jc,ipp))*udzq
!
!   22 second derivatives of vx
!
            dq32=(vx(kc,jmm,ic) &
                 -2.0*vx(kc,jc,ic) &
                 +vx(kc,jpp,ic))*udyq


          qcap(kc,jc,ic) =-(h31+h32+h33)+dq31+dq32+tempit
            
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end
!
