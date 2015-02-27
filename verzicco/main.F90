      program AFiD
      use mpih
      use param
      use local_arrays, only: vx,vy,vz,temp,pr
      use hdf5
      use decomp_2d
      use decomp_2d_fft
      use stat_arrays, only: nstatsamples
!$    use omp_lib
      implicit none
      integer :: errorcode, nthreads
      real    :: instCFL,dmax
      real    :: ti(2), tin(3), minwtdt
      real :: ts
      integer :: prow=0,pcol=0
      character(100) :: arg

!*******************************************************
!******* Read input file bou.in by all processes********
!*******************************************************
!
      call ReadInputFile

      if (command_argument_count().eq.2) then
        call get_command_argument(1,arg)
        read(arg,'(i)')prow
        call get_command_argument(2,arg)
        read(arg,'(i)')pcol
      endif

      call decomp_2d_init(nxm,nym,nzm,prow,pcol, &
     & (/ .false.,.true.,.true. /))

      ts=MPI_WTIME()
      tin(1) = MPI_WTIME()

      call MpiBarrier

      call HdfStart

      if (nrank.eq.master) ismaster = .true.

      if (ismaster) write(6,*) 'MPI tasks=', nproc

!$    if (ismaster) then 
!$OMP PARALLEL
!$OMP MASTER
!$        nthreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
!$        write(6,*) 'OMP threads=', nthreads
!$    end if

      if(ismaster) then
!m==========================================    
      call ResetLogs
!m====================================================
      write(6,112)ylen/alx3,zlen/alx3
  112 format(//,20x,'R A Y B E N ',//,10x, &
       '3D Cell with aspect-ratio:  D1/H = ',f5.2,' D2/H = ',f5.2)
      write(6,142) 
  142 format(//,8x,'Periodic lateral wall boundary condition')
      write(6,202) ray,pra
  202 format(/,5x,'Parameters: ',' Ra=',e10.3,' Pr= ',e10.3)
      if(variabletstep) then
         write(6,204) limitCFL
  204 format(/,5x,'Variable dt and fixed cfl= ', &
       e11.4,/ )            
      else 
         write(6,205) dtmax,limitCFL
  205 format(/,5x,'Fixed dt= ',e11.4,' and maximum cfl=', &
        e11.4,/ )            
      endif
      endif

      call InitTimeMarchScheme

      call InitVariables

      call CreateGrid

      call WriteGridInfo

      if (dumpslabs) call InitializeSlabDump

!m===================================                                                      
!m===================================
!m===================================
      if(ismaster) then
      write(6,754)nx,ny,nz                                              
  754 format(/,5x,'grid resolution: ',' nx= ',i5,' ny= ',i5, &
      ' nz= ',i5/)                       
      write(6,755) 1.d0/dx,1.d0/dy,1.d0/dz,dt,ntst                  
  755 format(/,2x,' dx=',e10.3,' dy=',e10.3,' dz=',e10.3,' dt=' &
      ,e10.3,' ntst=',i7,/)
      endif

!m===================================
!m===================================     
      
      time=0.d0
      if(statcal) nstatsamples = 0

      call InitPressureSolver
      call SetTempBCs

      if(readflow) then

        if(ismaster) write(6,*) 'Reading initial condition from file'

        call ReadFlowField

      else

        if(ismaster) write(6,*) 'Creating initial condition'

        ntime=0                                                         
        time=0.d0
        instCFL=0.d0
        
        call CreateInitialConditions

      endif                                                             

!EP   Update all relevant halos
      call update_halo(vx,1)
      call update_halo(vy,1)
      call update_halo(vz,1)
      call update_halo(temp,1)
      call update_halo(pr,1)

!EP   Check divergence. Should be reduced to machine precision after the first
!phcalc. Here it can still be high.

      call CheckDivergence(dmax)

      if(ismaster) write(6,*)' Initial maximum divergence: ',dmax

!EP   Write some values
      if(variabletstep) then
       if(ismaster) write(6,*)ntime,time,dt,dmax,tempm,tempmax,tempmin
      else
       if(ismaster) write(6,*)ntime,time,dt,dmax,tempm,tempmax,tempmin !RO Fix??
      end if
                                                                        
!  ********* starts the time dependent calculation ***
      errorcode = 0 !EP set errocode to 0 (OK)
      minwtdt = huge(0.0d0) !EP initialize minimum time step walltime
      do ntime=0,ntst                                           
        ti(1) = MPI_WTIME()

!EP   Determine timestep size
        call CalcMaxCFL(instCFL)

        if(variabletstep) then
          if(ntime.gt.1) then
            if(instCFL.lt.1.0d-8) then !EP prevent fp-overflow
              dt=dtmax
            else
              dt=limitCFL/instCFL
            endif
            if(dt.gt.dtmax) dt=dtmax
          else
            dt=dtmin
          endif
            if(dt.lt.dtmin) errorcode = 166
        else  
!RO    fixed time-step
          instCFL=instCFL*dt
          if(instCFL.gt.limitCFL) errorcode = 165
        endif

        call TimeMarcher

        if(mod(time,tout).lt.dt) then
         if(ismaster) then
          write(6,*) ' ---------------------------------------- '
          write(6,*) ' T = ',time,' NTIME = ',ntime,' DT = ',dt
         endif
        endif

        time=time+dt

        if(ntime.eq.1.or.mod(time,tout).lt.dt) then
          call GlobalQuantities
          if(vmax(1).gt.limitVel.and.vmax(2).gt.limitVel) errorcode = 266

            call CalcMaxCFL(instCFL)
            call CheckDivergence(dmax)
            call CalcPlateNu

            if(time.gt.tsta) then

             if (statcal)  call CalcStats
             if (dumpslabs) call SlabDumper
             if (disscal.and.statcal) call CalcDissipationNu

            endif

            if(.not.variabletstep) instCFL=instCFL*dt

            if(dmax.gt.resid) errorcode = 169

        endif

        if(time.gt.tmax) errorcode = 333

        ti(2) = MPI_WTIME()
        minwtdt = min(minwtdt,ti(2) - ti(1))
        if(mod(time,tout).lt.dt) then
          if(ismaster) then
          write(6,*) 'Maximum divergence = ', dmax
          write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),dmax,tempm,tempmax,tempmin
          write(6,*) 'Minimum Iteration Time = ', minwtdt, ' sec.'
          endif
          minwtdt = huge(0.0d0)
        endif

       if( (ti(2) - tin(1)) .gt. walltimemax) errorcode = 334


!EP   Conditional exits
      if(errorcode.ne.0) then

!EP    dt too small
        if(errorcode.eq.166) call QuitRoutine(tin,.false.,errorcode)

!EP   cfl too high    
        if(errorcode.eq.165) call QuitRoutine(tin,.false.,errorcode)
      
!EP   velocities diverged
        if(errorcode.eq.266) call QuitRoutine(tin,.false.,errorcode)
          
!EP   mass not conserved
        if(errorcode.eq.169) call QuitRoutine(tin,.false.,errorcode)

!EP   Physical time exceeded tmax, no error; normal quit
        if(errorcode.eq.333) call QuitRoutine(tin,.true.,errorcode)

!EP   walltime exceeded walltimemax, no error; normal quit
        if(errorcode.eq.334) call QuitRoutine(tin,.true.,errorcode)

       endif

      enddo !EP main loop

      call QuitRoutine(tin,1,errorcode)
      
      end                                                               

