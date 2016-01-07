! ---------------------------------------------------------------
! This routine factorizes the grid resolution in the horizontal 
! directions to make sure that efficient FFT's are used. The input
! is (In2) is "ny or ny" and the output (Out2) is the largest factor
! A loop in the main routine checks whether the largest factor is 
! smaller than 7. If not the simulation is aborted.
! ---------------------------------------------------------------
      subroutine  Factorize(In2,Out2)
      implicit none
      integer,intent(in)  :: In2  ! Horizontal grid resolution 
      integer,intent(out) :: Out2 ! Largest factor
      integer  :: Input,Divisor
      
      Input=In2
      ! Remove all factors of 2
      do 
        if (mod(Input,2) /= 0 .or. Input == 1) exit 
        Input = Input / 2       ! remove this factor 
      enddo

      ! consider odd factors
      Divisor = 3 
      do ! try 3, 5, 7, etc
      if (Divisor > Input) exit  ! if a factor is too large, exit and done
      do                         ! try this factor repeatedly
        if (mod(Input,Divisor) /= 0 .or. Input == 1) exit 
        Input = Input / Divisor   ! remove this factor 
      enddo
      Divisor = Divisor + 2   ! Next odd number
      enddo
      Out2=Divisor ! Largest factor. Should be small for effcient FFT

      return 
      end
