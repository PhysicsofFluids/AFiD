!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: interp.F90                                     !
!    CONTAINS: subroutine interp,gridnew,interptrilin     !
!                                                         ! 
!    PURPOSE: Trilinear interpolation to a new grid for   !
!     continuation files. Assumes alx3=1                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine interp(arrold,arrnew,nzo,nyo,nxo, &
     & istro3,stro3,intvar,xs2o,xe2o,xs3o,xe3o)
      use param
      use mpih
      use decomp_2d, only: xstart,xend
      implicit none
      integer :: istro3
      real :: stro3
      integer,intent(in) :: intvar,nyo,nxo,nzo
      integer,intent(in) :: xs2o,xe2o,xs3o,xe3o
 
      real,intent(out),dimension(1:nx,xstart(2):xend(2), &
     & xstart(3):xend(3)) :: arrnew
      real,dimension(0:nxo+1,xs2o-1:xe2o+1,xs3o-1:xe3o+1) :: arrold
      real,dimension(0:nzo+1) :: tcold,tmold
      real,dimension(1:nz) :: tcnew,tmnew
      real,dimension(0:nyo+1) :: rcold,rmold
      real,dimension(1:ny) :: rcnew,rmnew
      real,dimension(0:nxo+1) :: zzold,zmold
      real,dimension(1:nx) :: zznew,zmnew

      real, dimension(0:nzo+1) :: xold
      real, dimension(0:nyo+1) :: yold
      real, dimension(0:nxo+1) :: zold
      real, dimension(1:nz) :: xnew
      real, dimension(1:ny) :: ynew
      real, dimension(1:nx) :: znew

      real :: bn(6),an(8)

      real :: bix1,bix2,biy1,biy2,biz1,biz2,bix,biy,biz
      integer :: j,k,i,l
      real :: sfcf,sfff,sccf,scff,sfcc,sffc,sccc,scfc
      integer :: bici,bifi,bicj,bifj,bick,bifk

!EP   Create old grid
      call gridnew(nzo,zlen,0,1.0,tcold(1:nzo),tmold(1:nzo))
      call gridnew(nyo,ylen,0,1.0,rcold(1:nyo),rmold(1:nyo))
      call gridnew(nxo,1.0d0,istro3,stro3,zzold(1:nxo),zmold(1:nxo))

!EP   2nd order extrapolation of grid
      tcold(0) = 2*tcold(1)-tcold(2)
      tcold(nzo+1) = 2*tcold(nzo)-tcold(nzo-1)
      tmold(0) = 2*tmold(1)-tmold(2)
      tmold(nzo+1) = 2*tmold(nzo)-tmold(nzo-1)

      rcold(0) = 2*rcold(1)-rcold(2)
      rcold(nyo+1) = 2*rcold(nyo)-rcold(nyo-1)
      rmold(0) = 2*rmold(1)-rmold(2)
      rmold(nyo+1) = 2*rmold(nyo)-rmold(nyo-1)

      zzold(0) = 2*zzold(1)-zzold(2)
      zzold(nxo+1) = 2*zzold(nxo)-zzold(nxo-1)
      zmold(0) = 2*zmold(1)-zmold(2)
      zmold(nxo+1) = 2*zmold(nxo)-zmold(nxo-1)

!EP   Create new grid
      call gridnew(nz,zlen,0,1.0,tcnew,tmnew)
      call gridnew(ny,ylen,0,1.0,rcnew,rmnew)
      call gridnew(nx,1.0d0,istr3,str3,zznew,zmnew)
      
      select case (intvar)
          case (1)
             xold=tcold
             yold=rmold
             zold=zmold
             xnew=tcnew
             ynew=rmnew
             znew=zmnew
          case (2)
             xold=tmold
             yold=rcold
             zold=zmold
             xnew=tmnew
             ynew=rcnew
             znew=zmnew
          case (3)
             xold=tmold
             yold=rmold
             zold=zzold
             xnew=tmnew
             ynew=rmnew
             znew=zznew
          case (4)
             xold=tmold
             yold=rmold
             zold=zmold
             xnew=tmnew
             ynew=rmnew
             znew=zmnew
      end select

!EP   Extrapolate halo

      do i=xs3o,xe3o
        do j=xs2o,xe2o
          arrold(0    ,j,i)=2.0*arrold(1  ,j,i)-arrold(2,j,i)
          arrold(nxo+1,j,i)=2.0*arrold(nxo,j,i)-arrold(nxo-1,j,i)
        enddo
      enddo

      do j=xs2o,xe2o
        do k=1,nxo
          arrold(k,j,xs3o-1)=2.0*arrold(k,j,xs3o)-arrold(k,j,xs3o+1)
          arrold(k,j,xe3o+1)=2.0*arrold(k,j,xe3o)-arrold(k,j,xe3o-1)
        enddo
      enddo

      do i=xs3o,xe3o
        do k=1,nxo
          arrold(k,xs2o-1,i)=2.0*arrold(k,xs2o,i)-arrold(k,xs2o+1,i)
          arrold(k,xe2o+1,i)=2.0*arrold(k,xe2o,i)-arrold(k,xe2o-1,i)
        enddo
      enddo


!EP   INTERP
      do i=xstart(3),xend(3)
       do j=xstart(2),xend(2)
        do k=1,nx
!c   Find nearest grid value
      bix=xnew(i)
      biy=ynew(j)
      biz=znew(k)

      bifi=nzo-1
      bici=nzo
      do l=1,nzo
      if(xold(l).ge.bix) then
      bifi=l-1
      bici=l
      goto 10
      endif
      enddo
10    continue

      bifj=nyo-1
      bicj=nyo
      do l=1,nyo
      if(yold(l).ge.biy) then
      bifj=l-1
      bicj=l
      goto 20
      endif
      enddo
20    continue

      bifk=nxo-1
      bick=nxo
      do l=1,nxo
      if(zold(l).ge.biz) then
      bifk=l-1
      bick=l
      goto 30
      endif
      enddo
30    continue

!cc   Define
      bix1 = xold(bifi) 
      bix2 = xold(bici)
      biy1 = yold(bifj) 
      biy2 = yold(bicj)
      biz1 = zold(bifk) 
      biz2 = zold(bick)

!EP   Send data

       sfcf = arrold(bifk,bicj,bifi)
       sfff = arrold(bifk,bifj,bifi)
       sccf = arrold(bifk,bicj,bici)
       scff = arrold(bifk,bifj,bici)

       sfcc = arrold(bifk,bicj,bifi)
       sffc = arrold(bifk,bifj,bifi)
       sccc = arrold(bifk,bicj,bici)
       scfc = arrold(bifk,bifj,bici)

       an(1) = sfcf
       an(2) = sfff
       an(3) = sccf
       an(4) = scff
       an(5) = sfcc
       an(6) = sffc
       an(7) = sccc
       an(8) = scfc
       bn(1) = bix1
       bn(2) = bix2
       bn(3) = biy1
       bn(4) = biy2
       bn(5) = biz1
       bn(6) = biz2
       call interptrilin(bix,biy,biz,an,bn,arrnew(k,j,i))

        enddo
       enddo
      enddo

      end

      subroutine interptrilin(bix,biy,biz,an,bn,ans)
      implicit none
      
      real, intent(in) :: bix,biy,biz
      real,intent(in) :: bn(6),an(8)
      real, intent(out) :: ans

      real :: bix1,bix2,biy1,biy2,biz1,biz2
      real :: afifjck,acifjck,aficjck,acicjck
      real :: afifjfk,acifjfk,aficjfk,acicjfk,dxdydz
      

      aficjfk = an(1)
      afifjfk = an(2)
      acicjfk = an(3)
      acifjfk = an(4)
      aficjck = an(5)
      afifjck = an(6)
      acicjck = an(7)
      acifjck = an(8)
      
      bix1 = bn(1)
      bix2 = bn(2)
      biy1 = bn(3)
      biy2 = bn(4)
      biz1 = bn(5)
      biz2 = bn(6)

      dxdydz = 1.0/((bix2-bix1)*(biy2-biy1)*(biz2-biz1))

      ans = (aficjfk*(bix2-bix)*(biy-biy1)*(biz2-biz) &
     &      +afifjfk*(bix2-bix)*(biy2-biy)*(biz2-biz) &
     &      +acicjfk*(bix-bix1)*(biy-biy1)*(biz2-biz) &
     &      +acifjfk*(bix-bix1)*(biy2-biy)*(biz2-biz) &
     &      +aficjck*(bix2-bix)*(biy-biy1)*(biz-biz1) &
     &      +afifjck*(bix2-bix)*(biy2-biy)*(biz-biz1) &
     &      +acicjck*(bix-bix1)*(biy-biy1)*(biz-biz1) &
     &      +acifjck*(bix-bix1)*(biy2-biy)*(biz-biz1) &
     &      )*dxdydz

      end

      subroutine gridnew(n,rext,istr,str,rc,rm)
      implicit none
      double precision,intent(in) :: rext,str
      double precision,dimension(1:n),intent(out) :: rc,rm
      double precision,dimension(1:n) :: etaz
      double precision,dimension(1:n+400) :: etazm
      double precision :: x3,etain,delet,pi
      integer,intent(in) :: n,istr
      integer :: i,nxmo,nclip

      if (istr.eq.0) then
        do i=1,n
          x3=real(i-1)/real(n-1)
          rc(i)=rext*x3
        enddo
       endif

      if(istr.eq.6) then
      pi=2.0d0*asin(1.0d0)
      nclip = int(str)
      nxmo = n+nclip+nclip
      do i=1,nxmo
        etazm(i)=+cos(pi*(float(i)-0.5)/float(nxmo))
      end do
      do i=1,n
        etaz(i)=etazm(i+nclip)
      end do
      delet = etaz(1)-etaz(n)
      etain = etaz(1)
      do i=1,n
        etaz(i)=etaz(i)/(0.5*delet)
      end do
      rc(1) = 0.0
      do i=2,n-1
        rc(i) = rext*(1.-etaz(i))*0.5
      end do
      rc(n) = rext
      endif


      do i=1,n-1
        rm(i)=(rc(i)+rc(i+1))*0.5d0
      enddo
      rm(n) = 2*rc(n)-rm(n-1)

      return
      end
