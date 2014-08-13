      subroutine initia
      use param
      use local_arrays
      use decomp_2d, only: xstart,xend
      implicit none
      
      pr=0.d0
      rhs=0.d0
      ru1=0.d0
      ru2=0.d0
      ru3=0.d0
      ruro=0.d0
      dph=0.d0
      dphhalo=0.d0
      q1=0.d0
      q2=0.d0
      q3=0.d0
      dens=1.d0

      tc = 0.d0
      tm = 0.d0
      rc = 0.d0
      rm = 0.d0
      am3ssk = 0.d0
      ac3ssk = 0.d0
      ap3ssk = 0.d0
      zz = 0.d0
      zm = 0.d0
      g3rc = 0.d0
      g3rm = 0.d0

      return 
      end   
