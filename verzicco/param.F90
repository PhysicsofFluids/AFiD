!===========================================================
! Declaration of global variables
!***********************************************************      
      module param
        implicit none
!==========================================================			
!       read from input file bou.in
!==========================================================
        integer   :: nz, ny, nx
        integer   :: nsst, nwrit, nread, ntst, ireset
        real      :: walltimemax,tpin,tmax
        real      :: alx3,str3
        integer   :: istr3
        real      :: rext,rext2
        real      :: ray,pra,dt,resid,cflmax
        integer   :: inslws,inslwn
        integer   :: starea,tsta
        real      :: dtmin,dtmax,cfllim,vlim
        integer   :: nson,idtv
        real   :: tframe
!=================================================
!       end of input file
!=================================================
        real :: time
!******* Grid parameters**************************
        real :: dx2,dx3,dx1
        real :: dx2q,dx3q,dx1q
!        
        real, allocatable, dimension(:) :: tc,tm
        real, allocatable, dimension(:) :: rc,rm
        real, allocatable, dimension(:) :: zz,zm,g3rc,g3rm
!====================================================
!******* QUANTITIES FOR DERIVATIVES******************
        real, allocatable, dimension(:) :: udx3c,udx3m
!==========================================================
!******* Grid indices**************************************
        integer, allocatable, dimension(:) :: kmc,kpc,kmv,kpv
!===========================================================
!******* Metric coefficients *******************************
        real, allocatable, dimension(:) :: ap3ck,ac3ck,am3ck
        real, allocatable, dimension(:) :: ap3sk,ac3sk,am3sk
        real, allocatable, dimension(:) :: ap3ssk,ac3ssk,am3ssk   
!============================================================
!******* Variables for FFTW and Poisson solver****************
        real, allocatable, dimension(:) :: ak2,ap
        real, allocatable, dimension(:) :: ak1,ao
        real, allocatable, dimension(:) :: amphk,acphk,apphk
        
!===========================================================
!******* Other variables ***********************************
        integer  :: nzm, nym, nxm
        integer  :: iaxsy
        real :: rint
        real :: ren, pec
        real :: pi
        real :: al,ga,ro
        real :: beta
        real :: qqmax,qqtot
        real :: re
        real :: denmax,denmin,densm
        integer :: ntime
        integer, parameter:: ndv=3
        real, dimension(1:ndv) :: vmax
        real, dimension(1:3) :: gam,rom,alm
        real, allocatable, dimension(:,:) :: denbs,denbn
              
        logical :: dumpslabs,statcal,disscal
        logical :: readflow, readstats, ismaster=.false.
        logical :: resetlogs,variabletstep

      end module param
      
!************* End of param module******************************
!===============================================================
!******* 2D arrays, dynamically allocated by each process*******
      module local_arrays
      use param
        implicit none
        real,allocatable,dimension(:,:,:) :: q1,q2,q3
        real,allocatable,dimension(:,:,:) :: pr,dens,rhs
        real,allocatable,dimension(:,:,:) :: ru1,ru2,ru3,ruro
        real,allocatable,dimension(:,:,:) :: dph,qcap,dq,hro,dphhalo
      end module local_arrays

!===============================================================
      module stat_arrays
       implicit none
       real,allocatable, dimension(:) :: q1_me,q1_rms 
       real,allocatable, dimension(:) :: q2_me,q3_me,q2_rms,q3_rms 
       real,allocatable, dimension(:) :: dens_me,dens_rms 
       real, allocatable,dimension(:) :: disste,dissth,densq3_me
       integer :: timeint_cdsp
      end module stat_arrays
!=====================================================       
      module stat3_param
        implicit none
        integer :: kslab(1:9)
        real    :: zslab(1:9)
      end module stat3_param
!=====================================================       
      module mpih
        implicit none
        include 'mpif.h'
        integer :: ierr
        integer, parameter :: master=0
        integer :: MDP = MPI_DOUBLE_PRECISION
      end module mpih
!====================================================
      module fftw_params
!        use param, only: m2m,m2mh,m1m
        use iso_c_binding

        type, bind(C) :: fftw_iodim
           integer(C_INT) n, is, os
        end type fftw_iodim

        interface
          type(C_PTR) function fftw_plan_guru_dft(rank,dims, &
           howmany_rank,howmany_dims,in,out,sign,flags) &
           bind(C, name='fftw_plan_guru_dft')
           import
           integer(C_INT), value :: rank
           type(fftw_iodim), dimension(*), intent(in) :: dims
           integer(C_INT), value :: howmany_rank
           type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
           complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
           complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
           integer(C_INT), value :: sign
           integer(C_INT), value :: flags
         end function fftw_plan_guru_dft

           type(C_PTR) function fftw_plan_guru_dft_r2c(rank,dims, &
             howmany_rank,howmany_dims,in,out,flags) &
             bind(C, name='fftw_plan_guru_dft_r2c')
             import
             integer(C_INT), value :: rank
             type(fftw_iodim), dimension(*), intent(in) :: dims
             integer(C_INT), value :: howmany_rank
             type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
             real(C_DOUBLE), dimension(*), intent(out) :: in
             complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
             integer(C_INT), value :: flags
           end function fftw_plan_guru_dft_r2c
           
           type(C_PTR) function fftw_plan_guru_dft_c2r(rank,dims, &
             howmany_rank,howmany_dims,in,out,flags)  &
             bind(C, name='fftw_plan_guru_dft_c2r')
             import
             integer(C_INT), value :: rank
             type(fftw_iodim), dimension(*), intent(in) :: dims
             integer(C_INT), value :: howmany_rank
             type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
             complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
             real(C_DOUBLE), dimension(*), intent(out) :: out
             integer(C_INT), value :: flags
           end function fftw_plan_guru_dft_c2r


       end interface

        integer FFTW_PATIENT, FFTW_FORWARD, FFTW_BACKWARD,FFTW_ESTIMATE
        parameter (FFTW_PATIENT=32)   
        parameter (FFTW_ESTIMATE=64)   
        parameter (FFTW_FORWARD=-1)   
        parameter (FFTW_BACKWARD=1)   
        type(C_PTR) :: fwd_guruplan_y,bwd_guruplan_y 
        type(C_PTR) :: fwd_guruplan_z,bwd_guruplan_z
        logical :: planned=.false.

      end module fftw_params
