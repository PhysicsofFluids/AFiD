Code structure
--------------

The main loop with the integration scheme can be found in TimeMarcher.F90. It calls the following routines:

 * ExplicitTermsVX
 * ExplicitTermsVY
 * ExplicitTermsVZ
 * ExplicitTermsTemp
 * ImplicitAndUpdateVX
 * ImplicitAndUpdateVY
 * ImplicitAndUpdateVZ
 * ImplicitAndUpdateTemp
 * update_halo 2x
 * CalculateLocalDivergence
 * SolvePressureCorrection
 * update_halo
 * CorrectVelocity
 * CorrectPressure
 * update_halo 5x 

List of input and output files
------------------------------

Input files (.in)
 *  bou.in - Detailed information of the input variables is given in the next section
 *  stst3.in 

Output data files (.h5)
 * continua_vx.h5
 * continua_vy.h5
 * continua_vz.h5
 * continua_temp.h5
 * continua_master.h5
 * cordin_info.h5 

Output log files (.out)
 * axicor.out
 * fact3.out
 * nu_plate.out
 * nu_vol.out
 * rms_vel.out

Input variables in bou.in
-------------------------
```
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bou.in
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NXM             : Number of computational cells in the wall normal direction
NYM, NZM        : Number of computational cells in the homogenous/periodic directions
NSST            : if (NSST = 3) third order Runge Kutta integrator
                  if (NSST = 1) Adams Bashforth integrator
NREAD           : if (NREAD = 1) Start simulation by reading in continuation files
                  if (NREAD = 0) Start simulation without continuation files, initial
                                 conditions created in the code 

NTST            : Total number of time steps
WALLTIMEMAX     : Maximum wall time of the simulation
TOUT            : Simulation time interval for calculating statistics
IRESET          : if (IRESET = 1) Time in simulation resets to 0
                  if (IRESET = 0) Time in simulation unchanged

ALX3D           : Height of wall bounded direction (always set to 1.0)
ISTR3           : Type of clustering of computational cells in wall bounded direction
                  if (ISTR3 = 0) Uniform clustering
                  if (ISTR3 = 4) Hyperbolic tangent-type clustering
                  if (ISTR3 = 6) Clipped Chebychev-type clustering
STR3            : Clipping parameter for non-uniform clustering (ISTR3 > 0)

YLEN, ZLEN      : Length of homogenous/periodic directions 

RA              : Rayleigh number
PRA             : Prandtl number
DT              : Time step to start the simulation
RESID           : Maximum allowable residual for mass conservation
CFLMAX          : Maximum CFL number for simulation (set to 1.2)

STATON          : if (STATON = 1) Calculate statistics, else set to 0
BALANCEON       : if (BALANCEON = 1) Calculates nusselt number through global balance
                  equations relating dissipation and heat transport
TSTA            : Time to start the statistics
STAREAD         : Read in statistics from continuation files

INSLWS          : if (INSLWS = 1) No-slip boundary condition on lower wall
                  if (INSLWS = 0) Free-slip boundary condition on lower wall  
INSLWN          : if (INSLWN = 1) No-slip boundary condition on upper wall
                  if (INSLWN = 0) Free-slip boundary condition on upper wall

IDTV            : if (IDTV = 1) variable time stepping
                  if (IDTV = 0) fixed time stepping
DTMIN           : Minimum time step for variable time stepping
DTMAX           : Maximum time step for variable time stepping
VLIM            : Maximum allowable local velocity in simulation

SLABDUMP(STST3) : Output wall normal slabs at locations specified in stst3.in
                  If this option is active, you need to create a directory named 'stst3'
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
```
