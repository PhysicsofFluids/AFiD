        subroutine Checksum3DArray(qua,cksum)
        use param
        use mpih
        use decomp_2d, only: xstart, xend
        implicit none
        integer i,j,k

      real, intent(in), dimension(1:n3,xstart(2)-1:xend(2)+1, &
     & xstart(3)-1:xend(3)+1)::qua
        real, intent(out) :: cksum

        cksum = 0.0

        do i=xstart(3),xend(3)
         do j=xstart(2),xend(2)
          do k=1,n3m
           cksum=cksum+qua(k,j,i)
          enddo
         enddo
        enddo
 
        call MpiSumRealScalar(cksum)

        return 
        end

