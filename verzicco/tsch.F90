!************************************************************************
!
!           SUBROUTINE  TSCHEM
!
!   This subroutine manages the whole integration scheme.
!   The following equations are solved:          
!   
!    ~~     n
!   Q  -  Q                n         n       n-1   alp       2  ~~   n 
!  --------- = -alp*grad (P ) + gam*H + rho*H   + ----- nabla ( Q + Q )
!    d t                                          2 Re
!
!          i                           i               i
!   where H  are the nonlinear terms, P  the pressure Q  the velocities
!       ~~
!   and Q  the provisional non solenoidal velocity field.
!   The superscripts (~~, n, n-1) indicate the time step level. 
!                        n
!   The nonlinear terms H  are computed in the routines HDNL*, while
!   in the routines INVTR* are computed the remaining terms, updated
!   the non linear terms and inverted the equation to find the provisional
!   field at the new time step.
!       ~~
!   The Q  velocity field is projected onto a solenoidal field by a 
!   scalar Phi computed through the equation
!
!                         2            1          ~~
!                    nabla (Phi ) =  ------ div ( Q  )
!                                    alp dt
!
!   The right hand side of this equation is computed in the routine
!   DIVG, while the equation is solved in PHCALC.
!
!   In the routine UPDVP the solenoidal velocity field at the new time
!   step is then computed through
!
!                n+1  ~~
!               Q   = Q  - alt*dt grad (Phi)
!
!   Finally in the routine PRCALC is updated the pressure field
!
!                n+1   n        alp dt      2
!               P   = P + Phi - ------ nabla (Phi)
!                                2 Re
!
!   When the scalar field is computed (density, concentration,
!   temperature) the routines HDNLRO and INVTRRO are used. The same
!   strategy at the velocity field is used, except that the scalar
!   field does not need any correction.
!
!   All variables are located on a staggered grid with the velocities
!   on the faces of the computational cell and all the scalars at the
!   centre. This is important when terms belonging to different equations
!   are avaluated.
!
!   Further details of the scheme can be found in the paper
!   "A finite-difference scheme for three-dimensional incompressible
!    flows in cylindrical coordinates" by R. Verzicco and P. Orlandi
!    J. of Comp. Phys. 1996.
!
!
      subroutine tschem
      use param
      use local_arrays
      use mpih
      use decomp_2d
      implicit none
      integer :: ns
      integer :: j,k,i
#ifdef DEBUG
      integer :: kstartp
      real :: cksum2,cksum3,mck2,mck3,cksum1,mck1,mck4,cksum4
#endif
! 
!m      REAL timef !etime_,t(2)
!
!   TIME INTEGRATION : implicit viscous, 3rd order RK (Adams Bashfort)  
!                                                                       
      do ns=1,nsst                                                 
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)
#ifdef DEBUG        
        mck1=0.0d0
        mck2=0.0d0
        mck3=0.0d0
        do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
        do k=1,n3m
        mck1=mck1+q1(k,j,i)
        mck2=mck2+q2(k,j,i)
        mck3=mck3+q3(k,j,i)
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'i_last cksum'
        write(*,*) 'q1cksum= ',cksum1
        write(*,*) 'q2cksum= ',cksum2
        write(*,*) 'q3cksum= ',cksum3
        endif
        mck1=0.0d0
        mck2=0.0d0
        mck3=0.0d0
        do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
        do k=1,n3m
        mck1=mck1+q1(k,j,i)
        mck2=mck2+q2(k,j,i)
        mck3=mck3+q3(k,j,i)
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'j_last cksum'
        write(*,*) 'q1cksum= ',cksum1
        write(*,*) 'q2cksum= ',cksum2
        write(*,*) 'q3cksum= ',cksum3
        write(*,*) 'starting tsch'
        endif
#endif

#ifdef DEBUG 
        if(nrank.eq.0) then
        write(*,*) 'starting hdnl1'
        endif
#endif
        call hdnl1
#ifdef DEBUG 
        if(nrank.eq.0) then
        write(*,*) 'starting hdnl2'
        endif
#endif
        call hdnl2

#ifdef DEBUG        
        if(nrank.eq.0) then
        write(*,*) 'starting hdnl3'
        endif
#endif
        
        call hdnl3

#ifdef DEBUG        
        mck1=0.0d0
        mck2=0.0d0
        mck3=0.0d0
        do k=1,n3m
        do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
        mck1=mck1+dph(k,j,i)
        mck2=mck2+qcap(k,j,i)
        mck3=mck3+dq(k,j,i)
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'j_last cksum'
        write(*,*) 'dqcksum= ',cksum3
        write(*,*) 'dphcksum= ',cksum1
        write(*,*) 'qcapcksum= ',cksum2
        endif
        mck1=0.0d0
        mck2=0.0d0
        mck3=0.0d0
        do k=1,n3m
        do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
        mck1=mck1+dph(k,j,i)
        mck2=mck2+qcap(k,j,i)
        mck3=mck3+dq(k,j,i)
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'i_last cksum'
        write(*,*) 'dqcksum= ',cksum3
        write(*,*) 'dphcksum= ',cksum1
        write(*,*) 'qcapcksum= ',cksum2
        endif
#endif

#ifdef DEBUG        
        if(nrank.eq.0) then
        write(*,*) 'starting hdnlro'
        endif
#endif
        call hdnlro                         !!  "    "      "

#ifdef DEBUG        
        mck1=0.0d0
        do k=1,n3m
        do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
        mck1=mck1+hro(k,j,i)
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'j_last'
        write(*,*) 'hrosum= ',cksum1
        endif
        mck1=0.0d0
        do k=1,n3m
        do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
        mck1=mck1+hro(k,j,i)
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'i_last'
        write(*,*) 'hrosum= ',cksum1
        endif
#endif

#ifdef DEBUG 
        if(nrank.eq.0) then
        write(*,*) 'starting invtr1'
        endif
#endif
        call invtr1
#ifdef DEBUG 
        if(nrank.eq.0) then
        write(*,*) 'starting invtr2'
        endif
#endif
        call invtr2

#ifdef DEBUG        
        if(nrank.eq.0) then
        write(*,*) 'starting invtr3'
        endif
#endif
        
        call invtr3


      ! MAKE ONLY IP HALO
        call update_halo(q1,1)

      ! MAKE ONLY JP HALO

        call update_halo(q2,1)

#ifdef DEBUG        
        mck1=0.0d0
        mck2=0.0d0
        mck3=0.0d0
        mck4=0.0d0
        do k=1,n3m
        do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
        mck1=mck1+q1(k,j,i)
        mck2=mck2+q2(k,j,i)
        mck3=mck3+q3(k,j,i)
        mck4=mck4+dens(k,j,i)
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck4,cksum4,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'i_last_cksum'
        write(*,*) 'q1cksum= ',cksum1
        write(*,*) 'q2cksum= ',cksum2
        write(*,*) 'q3cksum= ',cksum3
        write(*,*) 'densckm= ',cksum4
        endif
        mck1=0.0d0
        mck2=0.0d0
        mck3=0.0d0
        mck4=0.0d0
        do k=1,n3m
        do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
        mck1=mck1+q1(k,j,i)
        mck2=mck2+q2(k,j,i)
        mck3=mck3+q3(k,j,i)
        mck4=mck4+dens(k,j,i)
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck4,cksum4,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'j_last_cksum'
        write(*,*) 'q1cksum= ',cksum1
        write(*,*) 'q2cksum= ',cksum2
        write(*,*) 'q3cksum= ',cksum3
        write(*,*) 'densckm= ',cksum4
        endif
#endif

        call divg 
#ifdef DEBUG        
        mck2=0.0d0
        do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
        do k=1,n3m
        mck2=mck2+abs(dph(k,j,i))
        enddo
        enddo
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'dphcksum= ',cksum2
        write(*,*) 'starting phcalcTP'
        endif
#endif
        
        call phcalc

!EP this copy can be avoided by changing transpose_x_to_y_real and
!transpose_y_to_x_real so these routines can handles arrays with
!halo. This copy is a defacto array temporary. Using inferred size
!arrays in the transpose calls results in 5 more of these. Time spent on
!this copy is 0.1% for 65^3 grid.
        do i=xstart(3),xend(3)
          do j=xstart(2),xend(2)
            do k=1,n3m
              dphhalo(k,j,i) = dph(k,j,i)
            enddo
          enddo
        enddo

        call update_halo(dphhalo,1)

#ifdef DEBUG        
        mck2=0.0d0
        do k=1,n3m
        do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
        mck2=mck2+abs(dph(k,j,i))
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'dphcksum= ',cksum2
        write(*,*) 'starting updvp'
        endif
#endif
        
        call updvp                 !! SOLENOIDAL VEL FIELD

        call update_halo(q1,1)
        call update_halo(q2,1)
        call update_halo(q3,1)

#ifdef DEBUG        
        mck1=0.0d0
        mck2=0.0d0
        mck3=0.0d0
        do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
        do k=1,n3m
        mck1=mck1+abs(q1(k,j,i))
        mck2=mck2+abs(q2(k,j,i))
        mck3=mck3+abs(q3(k,j,i))
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'q1cksum= ',cksum1
        write(*,*) 'q2cksum= ',cksum2
        write(*,*) 'q3cksum= ',cksum3
        endif
#endif
        
!       call update_both_ghosts(n1,n2,dens,kstart,kend)

        
        call prcalc                         !! PRESSURE FIELD

        call update_halo(pr,1)

#ifdef DEBUG        
        if(nrank.eq.0) then
        write(*,*) 'starting invtrro'
        endif
#endif
        
        call invtrro

#ifdef DEBUG        
        mck1=0.0d0
        do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
        do k=2,n3m
        mck1=mck1+abs(dens(k,j,i))
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'denscksum= ',cksum1
        endif
#endif

        call update_halo(dens,1)

#ifdef DEBUG        
        mck1=0.0d0
        do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
        do k=2,n3m
        mck1=mck1+abs(dens(k,j,i))
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'denscksum after up= ',cksum1
        endif
        mck1=0.0d0
        do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
        mck1=mck1+abs(dens(1,j,i))
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0, &
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'denscksum plat= ',cksum1
        endif
#endif

        enddo
!m================================       
!m================================
        if(mod(time,tpin).lt.dt) then
        if(nrank.eq.0) then
        write(6,*) ' ---------------------------------------- '
        write(6,*) ' T = ',time,' NTIME = ',ntime,' DT = ',dt
        endif
        endif
!m================================       
!m================================
      return                                                            
      end                                                               
!
