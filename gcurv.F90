      program gcurv
      use mpih
      use param
      use local_arrays, only: q2,q3,dens,pr,q1
      use hdf5
      use decomp_2d
      use decomp_2d_fft
#ifdef STATS
      use stat_arrays, only: timeint_cdsp
#endif
!$    use omp_lib
      implicit none
      integer :: ntstf, hdf_error, errorcode, nthreads
      real    :: cflm,dmax
      real    :: ti(2), tin(3)
      real :: ts

!*******************************************************
!******* Read input file bou.in by all processes********
!*******************************************************
!
      ts=MPI_WTIME()
      tin(1) = MPI_WTIME()

      call ReadInputFile

      call decomp_2d_init(n3m,n2m,n1m,0,0, &
     & (/ .false.,.true.,.true. /))

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

!m==========================================    
      call openfi
!m==========================================
!
      pi=2.d0*dasin(1.d0)                          
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

      call InitializeTimeMarchScheme

      call InitializeVariables

      call MakeGrid

      call h5open_f(hdf_error)


#ifdef STATS3
      call initstst3
#endif
!m===================================                                                      
#ifdef MOVIE
      call inimov
#endif
!m===================================
!m===================================
      if(nrank.eq.0) then
      write(6,754)n1,n2,n3                                              
  754 format(/,5x,'grid resolution: ',' n1= ',i5,' n2= ',i5, &
      ' n3= ',i5/)                       
      write(6,755) 1.d0/dx1,1.d0/dx2,1.d0/dx3,dt,ntst                  
  755 format(/,2x,' dx1=',e10.3,' dx2=',e10.3,' dx3=',e10.3,' dt=' &
      ,e10.3,' ntst=',i7,/)
      endif

!m===================================
!m===================================     
      
      time=0.d0

#ifdef STATS
!EP   Read or initialize stat arrays
      timeint_cdsp = 0
#endif

!EP   Initialize the pressure solver
      call phini
 

!EP   Set the temperature boundary conditions
      call densbo
!
!      create the initial conditions
!
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
711     format(3x,'check in cond : ','  ntstf =',i8,2x,'tpin =',f10.1//)
      endif

!EP   Write some values
      if(idtv.eq.1) then
        if(nrank.eq.0) then
          write(6,*)ntime,time,dt,dmax,densm,denmax,denmin
        endif
        else
        if(nrank.eq.0) then
          write(6,*)ntime,time,cflm,dmax, densm,denmax,denmin
        endif
        cflm=cflm*dt
      endif

      if(nrank.eq.0) then
        tin(2) = MPI_WTIME()
        write(6,*) 'Initialization Time = ', tin(2) -tin(1), ' sec.'
      endif
                                                                        
!  ********* starts the time dependent calculation ***
      errorcode = 0
      do ntime=0,ntstf                                           
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
          write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),dmax,densm,denmax,denmin
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
178         format(10x,'too large local residue for mass conservation: ',e12.5,' at ')
          endif
          call divgloc
          call quit(tin,0)
        endif

        if(errorcode.eq.333) then
!EP   Physical time exceeded tmax, no error; normal quit
          if(nrank.eq.0) then
            write(*,*) "time greater than tmax"
            write(*,*) "statistics and continuation updated"
          endif
          call quit(tin,1)
        endif

        if(errorcode.eq.334) then
!EP   walltime exceeded walltimemax, no error; normal quit
          if(nrank.eq.0) then
            write(*,*) "walltime greater than walltimemax"
            write(*,*) "statistics and continuation updated"
          endif
          call quit(tin,1)
        endif

        if(nrank.eq.0) then
          write(*,*) "unknown error"
        endif
        call quit(tin,0)

        endif

      enddo !EP main loop

      if(nrank.eq.0) then
        write(*,*) "Maximum number of timesteps reached"
      endif
      call quit(tin,1)
      
      stop                                                              
      end                                                               

