c***********************************************************************
!EP      Output the location(s) of excessive divergence
!EP      i.e. i,j,k where dqcap(i,j,k) > resid
c***********************************************************************

      subroutine divgloc
      use param
      use local_arrays, only: q2,q3,q1
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
            do kc=1,n3m
            kp=kc+1
              dqcap= (q1(kc,jc,ip)-q1(kc,jc,ic))*dx1
     &              +(q2(kc,jp,ic)-q2(kc,jc,ic))*dx2
     %              +(q3(kp,jc,ic)-q3(kc,jc,ic))*udx3m(kc)
              if (abs(dqcap).gt.resid) then
                write(*,*) ic,jc,kc,nrank
            write(*,*) "q1",(q1(kc,jc,ip)-q1(kc,jc,ic))*dx1
       write(*,*) "q2",(q2(kc,jp,ic)-q2(kc,jc,ic))*dx2
       write(*,*) "q3",(q3(kp,jc,ic)-q3(kc,jc,ic))*udx3m(kc)
                write(*,*) "q2m",ic,jc,kc,q2(kc,jc,ic)
                write(*,*) "q2p",ic,jp,kc,q2(kc,jp,ic)
              endif
      enddo
      enddo
      enddo
      
      return     
      end         
