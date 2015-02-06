      program AFiD
      use mpih
      use param
      use local_arrays, only: q2,q3,temp,pr,q1
      use hdf5
      use decomp_2d
      use decomp_2d_fft
      use stat_arrays, only: timeint_cdsp
!$    use omp_lib
      implicit none
      integer :: ntstf, errorcode, nthreads
      real    :: cflm,dmax
      real    :: ti(2), tin(3)
      real :: ts

!*******************************************************
!******* Read input file bou.in by all processes********
!*******************************************************
!
      call ReadInputFile

      call decomp_2d_init(nxm,nym,nzm,0,0, &
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
         write(6,204) cflmax
  204 format(/,5x,'Variable dt and fixed cfl= ', &
       e11.4,/ )            
      else 
         write(6,205) dtmax,cfllim
  205 format(/,5x,'Fixed dt= ',e11.4,' and maximum cfl=', &
        e11.4,/ )            
      endif
      endif

      call InitializeTimeMarchScheme

      call InitializeVariables

      call CreateGrid

      call WriteGridInfo

      if (dumpslabs) call InitializeSlabDump

!m===================================                                                      
!m===================================
!m===================================
      if(ismaster) then
      write(6,754)nz,ny,nx                                              
  754 format(/,5x,'grid resolution: ',' nx= ',i5,' ny= ',i5, &
      ' nz= ',i5/)                       
      write(6,755) 1.d0/dx1,1.d0/dx2,1.d0/dx3,dt,ntst                  
  755 format(/,2x,' dx1=',e10.3,' dx2=',e10.3,' dx3=',e10.3,' dt=' &
      ,e10.3,' ntst=',i7,/)
      endif

!m===================================
!m===================================     
      
      time=0.d0
      if(statcal) timeint_cdsp = 0

      call InitPressureSolver
      call SetTempBCs

      if(readflow) then

        if(ismaster) write(6,*) 'Reading initial condition from file'

        call ReadFlowField

      else

        if(ismaster) write(6,*) 'Creating initial condition'

        ntime=0                                                         
        time=0.d0
        cflm=0.d0
        
        call CreateInitialConditions

      endif                                                             

!EP   Update all relevant halos
      call update_halo(q1,1)
      call update_halo(q2,1)
      call update_halo(q3,1)
      call update_halo(temp,1)
      call update_halo(pr,1)

!EP   Check divergence. Should be reduced to machine precision after the first
!phcalc. Here it can still be high.

      call CheckDivergence(dmax)

      if(ismaster) write(6,*)' Initial maximum divergence: ',dmax

      ntstf=ntst                                                   

      if(ismaster) write(6,711) ntstf,tpin

711     format(3x,'check in cond : ','  ntstf =',i8,2x,'tpin =',f10.1//)

!EP   Write some values
      if(variabletstep) then
       if(ismaster) write(6,*)ntime,time,dt,dmax,tempm,tempmax,tempmin
      else
       if(ismaster) write(6,*)ntime,time,cflm,dmax, tempm,tempmax,tempmin
       cflm=cflm*dt
      endif

      if(ismaster) then
        tin(2) = MPI_WTIME()
        write(6,*) 'Initialization Time = ', tin(2) -tin(1), ' sec.'
      endif
                                                                        
!  ********* starts the time dependent calculation ***
      errorcode = 0
      do ntime=0,ntstf                                           
        ti(1) = MPI_WTIME()

!EP   Determine timestep size
        call CalcMaxCFL(cflm)

        if(variabletstep) then
          if(ntime.gt.1) then
            if(cflm.eq.0.0) then !EP prevent fp-overflow
              dt=dtmax
            else
              dt=cflmax/cflm
            endif
            if(dt.gt.dtmax) dt=dtmax
          else
            dt=dtmin
          endif
            if(dt.lt.dtmin) errorcode = 166
        else
          cflm=cflm*dt
          if(cflm.gt.cfllim) errorcode = 165
        endif

!EP   Integrate
        call TimeMarcher

        if(mod(time,tpin).lt.dt) then
         if(ismaster) then
          write(6,*) ' ---------------------------------------- '
          write(6,*) ' T = ',time,' NTIME = ',ntime,' DT = ',dt
         endif
        endif

        time=time+dt

        if(ntime.eq.1.or.mod(time,tpin).lt.dt) then
          call GlobalQuantities
          if(vmax(1).gt.vlim.and.vmax(2).gt.vlim) errorcode = 266

            call CalcMaxCFL(cflm)
            call CheckDivergence(dmax)
            call CalculatePlateNu

            if(time.gt.tsta) then

             if (statcal)  call CalcStats
             if (dumpslabs) call SlabDumper
             if (disscal.and.statcal) call CalcDissipationNu

            endif

            if(.not.variabletstep) cflm=cflm*dt

            if(dmax.gt.resid) errorcode = 169

        endif

        if(time.gt.tmax) errorcode = 333

        ti(2) = MPI_WTIME()
        if(mod(time,tpin).lt.dt) then
          if(ismaster) then
          write(6,*) 'Maximum divergence = ', dmax
          write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),dmax,tempm,tempmax,tempmin
          write(6,*) 'Iteration Time = ', ti(2) -ti(1), ' sec.'
          endif
        endif

       if( (ti(2) - tin(1)) .gt. walltimemax) errorcode = 334


!EP   Conditional exits
      if(errorcode.ne.0) then

!EP    dt too small
        if(errorcode.eq.166) call QuitRoutine(tin,0,errorcode)

!EP   cfl too high    
        if(errorcode.eq.165) call QuitRoutine(tin,0,errorcode)
      
!EP   velocities diverged
        if(errorcode.eq.266) call QuitRoutine(tin,0,errorcode)
          
!EP   mass not conserved
        if(errorcode.eq.169) call QuitRoutine(tin,0,errorcode)

!EP   Physical time exceeded tmax, no error; normal quit
        if(errorcode.eq.333) call QuitRoutine(tin,1,errorcode)

!EP   walltime exceeded walltimemax, no error; normal quit
        if(errorcode.eq.334) call QuitRoutine(tin,1,errorcode)

       endif

      enddo !EP main loop

      call QuitRoutine(tin,1,errorcode)
      
      end                                                               

