************************************************************************
c
c           SUBROUTINE  TSCHEM
c
c   This subroutine manages the whole integration scheme.
c   The following equations are solved:          
c   
c    ~~     n
c   Q  -  Q                n         n       n-1   alp       2  ~~   n 
c  --------- = -alp*grad (P ) + gam*H + rho*H   + ----- nabla ( Q + Q )
c    d t                                          2 Re
c
c          i                           i               i
c   where H  are the nonlinear terms, P  the pressure Q  the velocities
c       ~~
c   and Q  the provisional non solenoidal velocity field.
c   The superscripts (~~, n, n-1) indicate the time step level. 
c                        n
c   The nonlinear terms H  are computed in the routines HDNL*, while
c   in the routines INVTR* are computed the remaining terms, updated
c   the non linear terms and inverted the equation to find the provisional
c   field at the new time step.
c       ~~
c   The Q  velocity field is projected onto a solenoidal field by a 
c   scalar Phi computed through the equation
c
c                         2            1          ~~
c                    nabla (Phi ) =  ------ div ( Q  )
c                                    alp dt
c
c   The right hand side of this equation is computed in the routine
c   DIVG, while the equation is solved in PHCALC.
c
c   In the routine UPDVP the solenoidal velocity field at the new time
c   step is then computed through
c
c                n+1  ~~
c               Q   = Q  - alt*dt grad (Phi)
c
c   Finally in the routine PRCALC is updated the pressure field
c
c                n+1   n        alp dt      2
c               P   = P + Phi - ------ nabla (Phi)
c                                2 Re
c
c   When the scalar field is computed (density, concentration,
c   temperature) the routines HDNLRO and INVTRRO are used. The same
c   strategy at the velocity field is used, except that the scalar
c   field does not need any correction.
c
c   All variables are located on a staggered grid with the velocities
c   on the faces of the computational cell and all the scalars at the
c   centre. This is important when terms belonging to different equations
c   are avaluated.
c
c   Further details of the scheme can be found in the paper
c   "A finite-difference scheme for three-dimensional incompressible
c    flows in cylindrical coordinates" by R. Verzicco and P. Orlandi
c    J. of Comp. Phys. 1996.
c
c
      subroutine tschem
      use param
      use local_arrays
      use mpih
      use decomp_2d
      implicit none
      integer :: ns
#ifdef DEBUG
      integer :: kstartp
      real :: cksum2,cksum3,mck2,mck3,cksum1,mck1,mck4,cksum4
      integer :: j,k,i
#endif
      
      real :: cksum
c 
cm      REAL timef !etime_,t(2)
c
c   TIME INTEGRATION : implicit viscous, 3rd order RK (Adams Bashfort)  
c                                                                       
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
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0,
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
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0,
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
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0,
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
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0,
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
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0,
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
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0,
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
        call update_halo(q1,1,.TRUE.)

      ! MAKE ONLY JP HALO

        call update_halo(q2,1,.TRUE.)

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
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck4,cksum4,1,MDP,MPI_SUM,0,
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
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck4,cksum4,1,MDP,MPI_SUM,0,
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
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'dphcksum= ',cksum2
        write(*,*) 'starting phcalcTP'
        endif
#endif
        
        call phcalc
!       call update_both_ghosts(n1,n2+1,dph,kstart,kend)

        call update_halo(dph,1,.TRUE.)

#ifdef DEBUG        
        mck2=0.0d0
        do k=1,n3m
        do j=xstart(2),xend(2)
        do i=xstart(3),xend(3)
        mck2=mck2+abs(dph(k,j,i))
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'dphcksum= ',cksum2
        write(*,*) 'starting updvp'
        endif
#endif
        
        call updvp                 !! SOLENOIDAL VEL FIELD

        call update_halo(q1,1,.TRUE.)
        call update_halo(q2,1,.TRUE.)
        call update_halo(q3,1,.TRUE.)

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
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck3,cksum3,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(mck2,cksum2,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'q1cksum= ',cksum1
        write(*,*) 'q2cksum= ',cksum2
        write(*,*) 'q3cksum= ',cksum3
        endif
#endif
        
!       call update_both_ghosts(n1,n2,dens,kstart,kend)

        
        call prcalc                         !! PRESSURE FIELD

        call update_halo(pr,1,.TRUE.)

#ifdef DEBUG        
        if(nrank.eq.0) then
        write(*,*) 'starting invtrro'
        endif
#endif
        
        call invtrro

        call update_halo(dens,1,.TRUE.)

#ifdef DEBUG        
        mck1=0.0d0
        do i=xstart(3),xend(3)
        do j=xstart(2),xend(2)
        do k=2,n3m
        mck1=mck1+abs(dens(k,j,i))
        enddo
        enddo
        enddo
        call MPI_REDUCE(mck1,cksum1,1,MDP,MPI_SUM,0,
     &      MPI_COMM_WORLD,ierr)
        if(nrank.eq.0) then
        write(*,*) 'denscksum= ',cksum1
        endif
#endif

        enddo
cm================================       
cm================================
        if(mod(time,tpin).lt.dt) then
        if(nrank.eq.0) then
        write(6,*) ' ---------------------------------------- '
        write(6,*) ' T = ',time,' NTIME = ',ntime,' DT = ',dt
        endif
        endif
cm================================       
cm================================
      return                                                            
      end                                                               
c
