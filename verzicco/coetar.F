
c************************************************************************
c                                                                       *
c ****************************** subrout coetar  ********************** *
c                                                                       *
c    this subroutine calculates the coefficients for the                *
c    integration in the radial direction with non-uniform coor. trasf.  *
c                                                                       *
c************************************************************************
      subroutine coetar
      use param
      implicit none
      integer :: km,kc,kp
      real:: a33,a33m,a33p
c
c
c  ***********  coefficients for q3   for x3 differentation
c  c means centered that is at k location
c
      am3ck(1)=0.d0
      ap3ck(1)=0.d0
      ac3ck(1)=1.d0
      am3ck(n3)=0.d0
      ap3ck(n3)=0.d0
      ac3ck(n3)=1.d0
      do kc=2,n3m
      km=kc-1
      kp=kc+1
      a33=dx3q/g3rc(kc)
      a33p=1.d0/g3rm(kc)
      a33m=1.d0/g3rm(km)
      ap3ck(kc)=a33*a33p
      am3ck(kc)=a33*a33m
      ac3ck(kc)=-(ap3ck(kc)+am3ck(kc))
      enddo
c
c
c  **coefficients for q1, q2 ap3sk,am3sk,ac3sk
c   ap3ssk,am3ssk,ac3ssk, psc   for x3 differentation
c  s means staggered that is at k+1/2 location
c
      do kc=2,n3m-1
      kp=kc+1
      km=kc-1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      ap3sk(kc)=a33p
      am3sk(kc)=a33m
      ac3sk(kc)=-(ap3sk(kc)+am3sk(kc))
      enddo
c    
c    lower wall  bound. conditions  indicated by inslws
c    differemtiation of sec. derivative at 1+1/2
c    
      kc=1
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      ap3sk(kc)=a33p
      am3sk(kc)=0.d0
      ac3sk(kc)=-(a33p+inslws*a33m*2.d0)
c     ac3ssk(kc)=-(a33p)
c    
c    upper wall  bound. conditions  indicated by inslws
c    differentiation of sec. derivative at n3-1/2
c    
      kc=n3m
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      am3sk(kc)=a33m
      ap3sk(kc)=0.d0
      ac3sk(kc)=-(a33m+inslwn*a33p*2.d0)
c     ac3ssk(kc)=-(a33m)
C

      am3ssk(1)=0.d0
      ap3ssk(1)=0.d0
      ac3ssk(1)=1.d0

      do kc=2,n3m
       kp=kc+1
       km=kc-1
       a33=dx3q/g3rc(kc)
       a33p=1.d0/g3rm(kc)
       a33m=1.d0/g3rm(km)
       ap3ssk(kc)=a33*a33p
       am3ssk(kc)=a33*a33m
       ac3ssk(kc)=-(ap3ssk(kc)+am3ssk(kc))
      enddo

      am3ssk(n3)=0.d0
      ap3ssk(n3)=0.d0
      ac3ssk(n3)=1.d0

      return
      end
c

