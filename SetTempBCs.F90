!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: densbo.F90                                     !
!    CONTAINS: subroutine densbo                          !
!                                                         ! 
!    PURPOSE: Initialization routine. Calcuates the       !
!     temperature boundary conditions at the plates       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SetTempBCs
      use param
      implicit none
      integer :: j,i

      do i=1,n1m
       do j=1,n2m
        denbn(j,i)=0.d0
        denbs(j,i)=1.d0
       enddo
      enddo

      return
      end
!
