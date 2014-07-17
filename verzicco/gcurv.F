************************************************************************
      subroutine gcurv
      use mpih
      use param
      use local_arrays, only: q2,q3,dens,pr,q1
      use hdf5
      use decomp_2d
      use decomp_2d_fft
#ifdef STATS
      use stat_arrays, only: timeint_cdsp
#endif
      implicit none
      integer :: ntstf, hdf_error, errorcode
      real    :: cflm,dmax
      real    :: ti(2), tin(3)
      real :: ts
c
c     Code for the computation of three-dimensional incompressible flows    
c     in cylindrical polar coordinates.                                             
c                                                                       
c     This code solves flow fields bounded in the x3 (axial) and 
c     x2 (radial) directions. All boundaries can be no-slip or free-slip
c     by setting the appropriate indices in the input file.
c     The geometry  (a cylindrical can) is given in the subroutine cordi.
c     The discretization is uniform in the axial (3) and azimuthal (1)
c     directions, while can be non-uniform in the radial direction.
c
c     The equations for the following variables
c                                                                       
c      q1=v(theta)    q2=v(r)*r     q3=v(zeta)                          
c                                                                       
c     are discretized by finite-difference schemes.                    
c     The introduction of the variable q2 is necessary to avoid the
c     problem of the singularity of the equations at the axis of symmetry
c     (r=0).
c
c     All spatial derivatives are discretized by central second-order
c     accurate finite-difference schemes including the non linear terms.                                   
c
c     The nonlinear terms are treated explicitly, while the viscous terms
c     are computed implicitly. This would lead to the inversion of a
c     large banded matrix, that however is avoided by introducing
c     a factored scheme bringing to the solution of three tridiagonal
c     matrices for each velocity component (subroutine INVTR*).
c                              
c     In time a fractional-step procedure is used in the version of 
c     Nagi Mansour introducing the pressure in the first step.                         
c
c     The non-linear terms and the cross derivatives of the viscous terms
c     are discretized by explicit  Adams-Bashfort or 3rd order Runge-Kutta
c     method (A. Wray, personal communication).                      
c
c     The scalar quantity Phi, which projects the provisional velocity field
c     onto  a divergence free field, is solved by a direct method. 
c     For the axial and azimuthal directions modified wave numbers coupled 
c     with trigonometric expansions (FFTs) are used. The equation is then
c     solved by simply inverting a tridiagonal matrix for the radial direction.
c     No explicit boundary conditions are necessary for this Poisson equation.      
c                                                                       
c     Other details of the scheme are given in the introduction of the  
c     subroutine TSCHEM
c
c
      ts=MPI_WTIME()
      tin(1) = MPI_WTIME()
      
!EP   Initialize 
      call initia
c                                                                       
c     grid information                                                 
c                                                                       
      call meshes
      call indic                                                        
      call cordin

!RO   Make fftw plans
      call mkfftplans

      call h5open_f(hdf_error)
#ifdef STATS
      call initstst
#endif
#ifdef STATS3
      call initstst3
#endif
cm===================================                                                      
#ifdef MOVIE
      call inimov
#endif
cm===================================
cm===================================
      if(nrank.eq.0) then
      write(6,754)n1,n2,n3                                              
  754 format(/,5x,'grid resolution: ',' n1= ',i5,' n2= ',i5,
     % ' n3= ',i5/)                       
      write(6,755) 1.d0/dx1,1.d0/dx2,1.d0/dx3,dt,ntst                  
  755 format(/,2x,' dx1=',e10.3,' dx2=',e10.3,' dx3=',e10.3,' dt='
     & ,e10.3,' ntst=',i7,/)
      endif

cm===================================
cm===================================     
      
      time=0.d0

!EP   Initialize the pressure solver
      call phini
 

!EP   Set the temperature boundary conditions
      call densbo
c
c      create the initial conditions
c
      if(nread.eq.0) then
        if(nrank.eq.0) then
          write(6,*)' nread = 0: creating initial condition'
        endif
!EP   Set times to 0
        ntime=0                                                         
        time=0.d0
        cflm=0.d0
        
!EP   Create an initial condition
        call inqpr
      else
        if(nrank.eq.0) then
          write(6,*)' nread = 1: reading initial condition from file'
        endif
!EP   Read (and interpolate) continuation files
        call inirea

!EP   Increase the maximum simulation time by the end time in the
!continuation files
        tmax = tmax + time
      endif                                                             

!EP   Update all relevant halos
      call update_halo(q1,1)
      call update_halo(q2,1)
      call update_halo(q3,1)
      call update_halo(dens,1)
      call update_halo(pr,1)

!EP   Check divergence. Should be reduced to machine precision after the first
!phcalc. Here it can still be high.
      call divgck(dmax)
      if(nrank.eq.0) then
        write(6,*)' initial maximum divergence: ',dmax
      endif

      ntstf=ntst                                                   
      if(nrank.eq.0) then
        write(6,711) ntstf,tpin
711     format(3x,'check in cond : ', 
     &       '  ntstf =',i8,2x,'tpin =',f10.1//)
      endif

!EP   Write some values
      if(idtv.eq.1) then
        if(nrank.eq.0) then
          write(6,*)ntime,time,
     &    dt,dmax,densm,denmax,denmin
        endif
        else
        if(nrank.eq.0) then
          write(6,*)ntime,time,
     %    cflm,dmax, densm,denmax,denmin
        endif
        cflm=cflm*dt
      endif

!EP   Initialize metrics
      call coetar

      if(nrank.eq.0) then
        tin(2) = MPI_WTIME()
        write(6,*) 'Initialization Time = ', tin(2) -tin(1), ' sec.'
      endif
                                                                        
c  ********* starts the time dependent calculation ***
      do ntime=0,ntstf                                           
        errorcode = 0
        ti(1) = MPI_WTIME()

!EP   Determine timestep size
        call cfl(cflm)

        if(idtv.eq.1) then
          if(ntime.ne.1) then
            dt=cflmax/cflm
!EP   Restrict dt
            if(dt.gt.dtmax) dt=dtmax
          endif
            if(dt.lt.dtmin) errorcode = 166
        else
          cflm=cflm*dt
          if(cflm.gt.cfllim) errorcode = 165
        endif
        beta=dt/ren*0.5d0

!EP   Integrate
        call tschem
        time=time+dt


        if(ntime.eq.1.or.mod(time,tpin).lt.dt) then
          call globalquantities
          if(vmax(1).gt.vlim.and.vmax(2).gt.vlim) errorcode = 266
            call cfl(cflm)
            call divgck(dmax)
            call densmc
            if(time.gt.tsta) then
#ifdef STATS
            call stst
#endif
#ifdef STATS3
            call stst3
#endif
#ifdef BALANCE
            call balance
#endif
            endif
            if(idtv.eq.0) then
              cflm=cflm*dt
            endif
            if(dmax.gt.resid) errorcode = 169

        endif

#ifdef MOVIE
        if(mod(time,tframe).lt.dt) then
          call mkmov
        endif
#endif
         
        if(time.gt.tmax) errorcode = 333

        ti(2) = MPI_WTIME()
        if(mod(time,tpin).lt.dt) then
          if(nrank.eq.0) then
          write(6,*) 'Maximum divergence = ', dmax
          write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),
     &    dmax,densm,denmax,denmin
          write(6,*) 'Iteration Time = ', ti(2) -ti(1), ' sec.'
          endif
        endif

       if( (ti(2) - tin(1)) .gt. walltimemax) errorcode = 334


!EP   Conditional exits
      if(errorcode.ne.0) then

        if(errorcode.eq.166) then
!EP    dt too small
          if(nrank.eq.0) then
            write(6,168) dt 
168         format(10x,'dt too small, DT= ',e14.7)
          endif
          call quit(tin,0)
        endif

        if(errorcode.eq.165) then
!EP   cfl too high    
          if(nrank.eq.0) then
            write(6,164) 
164         format(10x,'cfl too large  ')
          endif
          call quit(tin,0)
        endif
      
        if(errorcode.eq.266) then
!EP   velocities diverged
          if(nrank.eq.0) then
            write(6,268)
268         format(10x,'velocities diverged')
          endif
          call quit(tin,0)
        endif
          
        if(errorcode.eq.169) then
!EP   mass not conserved
          if(nrank.eq.0) then
            write(6,178) dmax                                 
178         format(10x,'too large local residue for mass conservation
     & : ',e12.5,' at ')
          endif
          call divgloc
          call quit(tin,0)
        endif

        if(errorcode.eq.333) then
!EP   Physical time exceeded tmax, no error; normal quit
          write(*,*) "time greater than tmax"
          write(*,*) "statistics and continuation updated"
          call quit(tin,1)
        endif

        if(errorcode.eq.333) then
!EP   walltime exceeded walltimemax, no error; normal quit
          write(*,*) "walltime greater than walltimemax"
          write(*,*) "statistics and continuation updated"
          call quit(tin,1)
        endif

        write(*,*) "unknown error"
        call quit(tin,0)

        endif

      enddo !EP main loop

      write(*,*) "Maximum number of timesteps reached"
      call quit(tin,1)
      
      return                                                            
      end                                                               

      subroutine quit(tin,cond)
      use hdf5
      use mpih
      use decomp_2d, only: nrank
      use decomp_2d_fft
      implicit none
      integer, intent(in) :: cond
      integer :: hdf_error
      real :: tin(3)

      tin(3) = MPI_WTIME()
      if(nrank.eq.0) then
          write(6,*) 'Total Iteration Time = ',tin(3) -tin(2),' sec.'
      endif

      call closefi
      call mem_dealloc
      call h5close_f(hdf_error)
      call decomp_2d_fft_finalize

      if(cond.eq.1) then
#ifdef STATS
        call ststwr
#endif
        call continua
      endif


      stop

      end subroutine quit
