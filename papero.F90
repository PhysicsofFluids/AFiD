!     this code is made for simulating two-dimensional flows in 
!     cartesian coordinates.
!     boundary condition are no-slip on the horizontal plates and either
!     periodic, no-slip or free-slip on the lateral wall
!                                                                       
!      navier-stokes equations are solved by a fractional step method   
!     ( Kim and Moin ) with the pressure in the first step.              
!            
!     The time advancement of the solution is obtained by a
!     Runge-Kutta 3rd order low storage scheme (Wray) or a 2nd   
!     order Adams-Bashfort scheme.                                      
!
!     The Poisson  equation for the pressure is solved directly                
!     introducing FFT in the azimutal and vertical direction.
!     Because of the presence of the walls in the vertical direction
!     the cosFFT in that direction is required
!                                                                       
!                 Roberto Verzicco and Paolo Orlandi                       
!                 dipartimento di meccanica ed aeronautica              
!                 universita' la sapienza di roma                 
!                                                                       
!     All variables are calculated in a staggered grid:
!   
!        q1=vx, q2= vy, q3= vz       
!
!        dq2,dq3 : velocity correction                              
!
!        qcap :divergence of the  non free divergent velocity field    
!
!        non linear terms:
!
!        ru2, ru3, ruro : old step      
!        h2,h3, hro : new step           
!        dens : density
!        pr : pressure                                                  
!        dph : pressure correction                                      
!       pressure solver is in the package press.f                      
!       non linear terms are calculated in hdnl routines                
!       the invertions of momentum equations is performed in invtr      
!       routines                                             
! 
!======================================================================
!     Cartesian version and hybrid parallelization by Erwin P. van der Poel, University of Twente
!======================================================================                                                                     
      program papero
      use mpih
      use decomp_2d
      use param
!$    use omp_lib

      implicit none
      character(len=4) :: dummy
      integer :: tfini,tfin,n,ns,nthreads

!*******************************************************
!******* Read input file bou.in by all processes********
!*******************************************************
!
      open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) n1,n2,n3,nsst,nwrit,nread
        read(15,301) dummy
        read(15,*) ntst,walltimemax,tpin,tmax,ireset
        read(15,301) dummy
        read(15,*) alx3,istr3,str3
        read(15,301) dummy
        read(15,*) rext,rext2
        read(15,301) dummy
        read(15,*) ray,pra,dt,resid,cflmax
        read(15,301) dummy
        read(15,*) tsta,starea
        read(15,301) dummy
        read(15,*) inslws,inslwn
        read(15,301) dummy       
        read(15,*) idtv,dtmin,dtmax,cfllim,vlim
        read(15,301) dummy       
        read(15,*) tframe
301     format(a4)                
      close(15)

!EP INIT
      n1m=n1-1                                                          
      n2m=n2-1                                                          
      n3m=n3-1

      call decomp_2d_init(n3m,n2m,n1m,0,0, &
     & (/ .false.,.true.,.true. /))

!EP   DEBUG:
!     call decomp_2d_init(n3m,n2m,n1m,2,1,
!    & (/ .false.,.true.,.true. /))

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (nrank .eq. 0)write(6,*) 'MPI tasks=', nproc

!$    if (nrank .eq. 0) then 
!$OMP PARALLEL
!$OMP MASTER
!$        nthreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
!$        write(6,*) 'OMP threads=', nthreads
!$    end if

!m============================================
!m************ End of input file**************
!m============================================
      if( n1>m1 ) then
      if(nrank.eq.0) then
          write(6,*) 'Error: n1 must be = m1'
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
      endif
      
      if( n2>m2 ) then
      if(nrank.eq.0) then
          write(6,*) 'Error: n2 must be = m2'
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
      endif
      
      if( n3>m3 ) then
      if(nrank.eq.0) then
          write(6,*) 'Error: n3 must be = m3'
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr )
      endif
      endif
!m============================================
      rint = 0.d0
!
!     DEFINITIONS FOR THE NATURAL CONVECTION
!
      ren = dsqrt(ray/pra)
      pec = dsqrt(pra*ray)
!                                                                       
!
!m==========================================    
      call openfi
!m==========================================
!
      pi=2.d0*dasin(1.d0)                          
      tfini=dt*ntst                                                     
      tfin=tfini                                                        
!m======================                                                          
!
!m====================================================
      if(nrank.eq.0) then
!m====================================================                                                                             
      write(6,112)rext/alx3
  112 format(//,20x,'R A Y B E N ',//,10x, &
       '2D Cell with aspect-ratio:  D/H = ',f5.2)
      write(6,142) 
  142 format(//,8x,'Periodic lateral wall boundary condition')
      write(6,202) ray,pra
  202 format(/,5x,'Parameters: ',' Ra=',e10.3,' Pr= ',e10.3)
      if(idtv.eq.1) then
         write(6,204) cflmax
  204 format(/,5x,'Variable dt and fixed cfl= ', &
       e11.4,/ )            
      else 
         write(6,205) dtmax,cfllim
  205 format(/,5x,'Fixed dt= ',e11.4,' and maximum cfl=', &
        e11.4,/ )            
      endif
!m====================================================    
      endif
!m====================================================      
!                                                                       
!     assign coefficients for time marching schemes                     
!
      if(nsst.gt.1) then   
        gam(1)=8.d0/15.d0
        gam(2)=5.d0/12.d0
        gam(3)=3.d0/4.d0
        rom(1)=0.d0
        rom(2)=-17.d0/60.d0
        rom(3)=-5.d0/12.d0
!m======================================================
      if(nrank.eq.0) then
        write(6,100) (gam(n),n=1,nsst),(rom(n),n=1,nsst)
  100   format(/,5x,'The time scheme is a III order Runge-Kutta' &
        ,4x,'gam= ',3f8.3,4x,'ro= ',3f8.3)
      endif
!m======================================================
      else                                                              
        gam(1)=1.5d0
        gam(2)=0.d0
        gam(3)=0.d0
        rom(1)=-0.5d0
        rom(2)=0.d0
        rom(3)=0.d0
!m======================================================
      if(nrank.eq.0) then
        write(6,110) gam(1),rom(1)
  110   format(/,5x,'The time scheme is the Adams-Bashfort',4x, &
         'gam= ',f8.3,4x,'ro= ',f8.3)
      endif
     
!m======================================================                                 
      endif                                                             
      do 10 ns=1,nsst
        alm(ns)=(gam(ns)+rom(ns))
   10 continue

!m======================================================
      
      call mpi_workdistribution

      call mem_alloc

!m======================================================

!
!     the solution of the problem starts
!
!m======================================================


      call gcurv
      
      stop                                                              
      end                                                               
