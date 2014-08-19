
!***********************************************************************
      subroutine meshes
      use param
      implicit none
      dx1=rext/real(n1m)
      dx2=rext2/real(n2m)
      dx3=alx3/real(n3m)
      dx1=1.0/dx1
      dx2=1.0/dx2
      dx3=1.0/dx3
      dx1q=dx1*dx1                                                      
      dx2q=dx2*dx2                                                      
      dx3q=dx3*dx3                                                      
      return                                                            
      end                                                               
