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
      integer :: ntstf, l, hdf_error
      real    :: cflm,dmax
      real    :: ti(2), tin(3)
      real :: ts,te
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
c     timings
c                                                                       
c
c
      ts=MPI_WTIME()
      tin(1) = MPI_WTIME()
      
c
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
#ifdef STATS2
      call initstst2
#endif
C**************************
cm       imovie = 0
cm       tframe = 0.05
C**************************
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
c                                                                       
c   read or create initial fields                                       
c                                                                       
#ifdef DEBUG
      write(*,*) 'starting phini'
#endif
      call phini
c
        do 22 l=1,ndv                                                   
           vmax(l)=0.d0
   22   continue                                                        
c
#ifdef DEBUG
      write(*,*) 'starting densbo'
#endif
      call densbo
c
c      create the initial conditions
c
      if(nread.eq.0) then
cm================================
        if(nrank.eq.0) then
        write(6,*)' nread=0 ---> cond. iniz. create nel code'     
      endif
cm================================            
        ntime=0                                                         
        time=0.d0
        cflm=0.d0
        
#ifdef DEBUG                                                               
        write(*,*) 'starting inqpr'
#endif
        call inqpr

        call update_halo(q1,1,.TRUE.)
        call update_halo(q2,1,.TRUE.)
        call update_halo(q3,1,.TRUE.)
        call update_halo(dens,1,.TRUE.)
        call update_halo(pr,1,.TRUE.)

        
cm======================================================================  
#ifdef DEBUG                                                               
        write(*,*) 'starting divgck' 
#endif
        call divgck(dmax)
        if(nrank.eq.0) then
        write(6,*)' initial divg dmax  ',dmax
      endif
cm================================ 
cm
      else
cm
cm================================
        if(nrank.eq.0) then
       write(6,*)' nread=1 ---> cond. iniz. lette da file'      
      endif
cm================================ 
cm      
#ifdef DEBUG                                                               
        write(*,*) 'starting inirea' 
#endif
       call inirea
        call update_halo(q1,1,.TRUE.)
        call update_halo(q2,1,.TRUE.)
        call update_halo(q3,1,.TRUE.)
        call update_halo(dens,1,.TRUE.)
        call update_halo(pr,1,.TRUE.)
      tmax = tmax + time
cm
cm=====================================================
       
#ifdef DEBUG                                                               
        write(*,*) 'starting divgck' 
#endif
       call divgck(dmax)
        if(nrank.eq.0) then
       write(6,*)' initial divg dmax  ',dmax
      endif
cm==============================               
      endif                                                             
c
      ntstf=ntst                                                   
cm================================
        if(nrank.eq.0) then
      write(6,711) tprint,ntstf,tpin
711   format(3x,'check in cond : tprint =',f10.1, 
     &       '  ntstf =',i8,2x,'tpin =',f10.1//)
      endif
cm================================ 
#ifdef DEBUG                                                               
        write(*,*) 'starting densmc' 
#endif
       call densmc
            if(idtv.eq.1) then
            
cm================================
        if(nrank.eq.0) then
      write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),
     & dt,dmax,densm,denmax,denmin
      endif
cm================================  
            else
cm================================
        if(nrank.eq.0) then
      write(6,*)ntime,time,vmax(1),vmax(2),vmax(3),
     % cflm,dmax, densm,denmax,denmin
      endif
cm================================     
     
        cflm=cflm*dt
            endif
c
#ifdef DEBUG                                                               
        write(*,*) 'starting coetar' 
#endif
         call coetar
c      
c
c  ********* starts the time dependent calculation ***

        if(nrank.eq.0) then
      tin(2) = MPI_WTIME()
      
       write(6,*) 'Initialization Time = ', tin(2) -tin(1), ' sec.'
       write(6,*) 'Max D',dmax
      endif
c                                                                       
      do 350 ntime=0,ntstf                                           
c
c     the calculation stops if the velocities are diverging for numerical
c     stability conditions (courant number restrictions)                
c
          ti(1) = MPI_WTIME()
c
#ifdef DEBUG                                                               
        if(nrank.eq.0) then
        write(*,*) 'starting cfl' 
        endif
#endif
            call cfl(cflm)
#ifdef DEBUG                                                               
        write(*,*) 'finished cfl' ,nrank,cflm
#endif
            
            
            if(idtv.eq.1) then
              if(ntime.ne.1) then
                 dt=cflmax/cflm
                 if(dt.gt.dtmax) dt=dtmax
              endif
              if(dt.lt. 0.00000001d0) go to 166
            else
              cflm=cflm*dt
              if(cflm.gt.cfllim) go to 165
            endif
            beta=dt/ren*0.5d0
#ifdef DEBUG                                                               
        write(*,*) 'starting tschem'
#endif
            call tschem
            time=time+dt
            if(ntime.eq.1) go to 306 
            if(mod(time,tpin).lt.dt) go to 306  
            go to 305
  306 continue                                                          
            call vmaxv
            if(vmax(1).gt.1000.d0.and.vmax(2).gt.1000.d0) go to 266
            call cfl(cflm)
            call divgck(dmax)
            call densmc
            if(time.gt.tsta) then
#ifdef STATS
            call stst
#endif
#ifdef STATS2
            call stst2
#endif
#ifdef BALANCE
            call balance
#endif
            endif
            if(idtv.eq.1) then
            else
              cflm=cflm*dt
            endif
            if(dmax.gt.resid) go to 169
c                                                                       
  305 continue                                                          

#ifdef MOVIE
       if(mod(time,tframe).lt.dt) then
        call mkmov
       endif
#endif
          
                 if(time.gt.tmax) go to 333
         ti(2) = MPI_WTIME()
              if(mod(time,tpin).lt.dt) then
              if(nrank.eq.0) then
               write(6,*) 'Iteration Time = ', ti(2) -ti(1), ' sec.'
          write(6,*) 'dmax',dmax
              endif
              endif
cm===============================================
           if( (ti(2) -tin(1)) .gt. 85500.) then
#ifdef STATS
        call ststwr
#endif
       call continua

         if(nrank.eq.0) then
            write(6,*) 'Restart file updated at exit'
            write(6,*) 'exit time=',ti(2) -tin(1),'sec'
         endif
         call MPI_FINALIZE(ierr)

         stop ! Exit after time limit

        endif
c================================================           
c
  350 continue
 
333   continue
  
        tin(3) = MPI_WTIME()
      if(nrank.eq.0) then
          write(6,*) 'Total Iteration Time = ',tin(3) -tin(2),' sec.'
      endif
#ifdef STATS
      call ststwr
#endif
      call continua
      go to 167                                                         
  165 continue 

cm======================================
      if(nrank.eq.0) then
      write(6,164)                                                      
  164 format(10x,'cfl too large  ')
      endif
cm======================================                                  
      go to 167                                                         
  166 continue
cm======================================
      if(nrank.eq.0) then
      write(6,168) dt 
  168 format(10x,'dt too small, DT= ',e14.7)
      endif
cm======================================   
      go to 167                                                         
  266 continue
cm======================================
      if(nrank.eq.0) then
      write(6,268)                                                      
  268 format(10x,'velocities diverged')
      endif
cm======================================                                 
      go to 167                                                         
  169 continue
  
cm======================================
      if(nrank.eq.0) then
      write(6,178) dmax                                 
  178 format(10x,'too large local residue for mass conservation : '     
     1       ,e12.5,' at ')
     
      endif
      call divgloc
cm======================================

  167 continue                                                          

      te=MPI_WTIME()
      
      if(nrank.eq.0) then
        open(27,file="Total_time.out")
        write(27,*)"Total simulation time in sec.: ", te-ts
        close(27)
      endif

      call mem_dealloc

      call h5close_f(hdf_error)
c
      call decomp_2d_fft_finalize
      
!     call MPI_FINALIZE(ierr)

      
      return                                                            
      end                                                               
c
