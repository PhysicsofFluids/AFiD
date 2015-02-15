!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: LocateLargeDivergence.F90                      !
!    CONTAINS: subroutine LocateLargeDivergence           !
!                                                         ! 
!    PURPOSE: Debugging routine. Output the location(s)   !
!     of excessive divergence.                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine LocateLargeDivergence
      use param
      use local_arrays, only: vy,vx,vz
      use mpih
      use decomp_2d, only: xstart,xend,nrank
      implicit none
      integer :: jc,kc,kp,jp,ic,ip
      real    :: dqcap
        
      if(nrank.eq.0) write(*,*) "I   J   K   RANK"
      do ic=xstart(3),xend(3)
        ip=ic+1
        do jc=xstart(2),xend(2)
          jp=jc+1
            do kc=1,nxm
            kp=kc+1
              dqcap= (vz(kc,jc,ip)-vz(kc,jc,ic))*dz &
     &              +(vy(kc,jp,ic)-vy(kc,jc,ic))*dy &
     &              +(vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)
              if (abs(dqcap).gt.resid) then
                write(*,*) ic,jc,kc,nrank
            write(*,*) "vz",(vz(kc,jc,ip)-vz(kc,jc,ic))*dz
       write(*,*) "vy",(vy(kc,jp,ic)-vy(kc,jc,ic))*dy
       write(*,*) "vx",(vx(kp,jc,ic)-vx(kc,jc,ic))*udx3m(kc)
                write(*,*) "vym",ic,jc,kc,vy(kc,jc,ic)
                write(*,*) "vyp",ic,jp,kc,vy(kc,jp,ic)
             endif
      enddo
      enddo
      enddo
      
      return     
      end         
