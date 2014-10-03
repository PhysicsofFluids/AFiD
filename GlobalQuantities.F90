!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: GlobalQuantities.F90                           !
!    CONTAINS: subroutine GlobalQuantities                !
!                                                         ! 
!    PURPOSE: Calculate maximum velocity and temperature. !
!     volume averaged Nusselt number and Reynolds number  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine GlobalQuantities
      use param
      use local_arrays, only: q2,q3,q1,temp
      use decomp_2d, only: xstart,xend
      use mpih
      implicit none
      integer :: jc,kc,kp,ic
      real :: anusin,vol,q3cen,fac2,tempcen
      real :: q2_rms_vol,q1_rms_vol
      real :: q3_rms_vol,q1q2q3_rms_vol,rradpr

!EP   Initialize
      vmax(1)=-huge(0.0d0)
      vmax(2)=-huge(0.0d0)
      vmax(3)=-huge(0.0d0)
      tempmax=-huge(0.0d0)
      tempmin=huge(0.0d0)
      tempm=0.0d0
      anusin=0.d0 
      q1_rms_vol = 0.0d0
      q2_rms_vol = 0.0d0
      q3_rms_vol = 0.0d0
      q1q2q3_rms_vol = 0.0d0
      vmax = 0.0d0
      vol = 1.d0/(alx3*dx3*real(nzm)*real(nym))


!EP   Loop over volume
      do ic=xstart(3),xend(3)
        do jc=xstart(2),xend(2)
          do kc=1,nxm
          kp = kc + 1
          fac2 = g3rc(kc)
          vmax(1) = max(vmax(1),abs(q1(kc,jc,ic)))
          vmax(2) = max(vmax(2),abs(q2(kc,jc,ic)))
          vmax(3) = max(vmax(3),abs(q3(kc,jc,ic)))
          tempmax = max(tempmax,temp(kc,jc,ic))
          tempmin = min(tempmin,temp(kc,jc,ic))
          q3cen = (q3(kc,jc,ic)+q3(kp,jc,ic))*0.5d0
          tempcen = (temp(kc,jc,ic)+temp(kp,jc,ic))*0.5d0
          anusin=anusin+tempcen*q3cen*fac2
          tempm=tempm+tempcen*fac2
          q1_rms_vol = q1_rms_vol + fac2*q1(kc,jc,ic)**2
          q2_rms_vol = q2_rms_vol + fac2*q2(kc,jc,ic)**2
          q3_rms_vol = q3_rms_vol + fac2*q3(kc,jc,ic)**2
          q1q2q3_rms_vol = q1q2q3_rms_vol + fac2* &
     &    (q1(kc,jc,ic)**2+q2(kc,jc,ic)**2+q3(kc,jc,ic)**2)  
          enddo
        enddo
      enddo

!EP   Reduce

      call MpiSumRealScalar(tempm)
      call MpiSumRealScalar(anusin)
      call MpiSumRealScalar(q1_rms_vol)
      call MpiSumRealScalar(q2_rms_vol)
      call MpiSumRealScalar(q3_rms_vol)
      call MpiSumRealScalar(q1q2q3_rms_vol)
      call MpiMinRealScalar(tempmin)
      call MpiMaxRealScalar(tempmax)
      call MpiMaxRealScalar(vmax(1))
      call MpiMaxRealScalar(vmax(2))
      call MpiMaxRealScalar(vmax(3))
       
!EP   Write logs
      if(ismaster) then

      anusin=1.d0 + dsqrt(pra*ray)*anusin*vol

      open(95,file='nu_vol.out',status='unknown',access='sequential', &
        position='append')
      write(95,*) time, anusin
      close(95)

      rradpr=dsqrt(ray/pra)
      q1_rms_vol=dsqrt(q1_rms_vol*vol)*rradpr
      q2_rms_vol=dsqrt(q2_rms_vol*vol)*rradpr
      q3_rms_vol=dsqrt(q3_rms_vol*vol)*rradpr
      q1q2q3_rms_vol=dsqrt(q1q2q3_rms_vol*vol)*rradpr

       open(94,file='rms_vel.out',status='unknown',position='append', &
        access='sequential')
       write(94,*) time,q1_rms_vol,q2_rms_vol,q3_rms_vol, &
     & q1q2q3_rms_vol
       close(94)

      endif

      return   
      end     
