c===========================================================
c Declaration of global variables
c***********************************************************      
      module param
        implicit none
        integer :: m2,m3,m1
        integer :: m2m,m3m,m2mh,m1m
c       Grid size        
       parameter (m1=33,m2=33,m3=65)
        parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1,m2mh=m2m/2+1)
c==========================================================			
c       read from input file bou.in
c==========================================================
        integer   :: n2, n3,n1
        integer   :: nsst, nwrit, nread, ntst, ireset
        real      :: tprint,tpin,tmax
        real      :: alx3,str3
        integer   :: istr3
        real      :: rext,rext2
        real      :: ray,pra,dt,resid,cflmax
        integer   :: inslws,inslwn
        integer   :: starea,tsta
        real      :: dtmax,cfllim 
        integer   :: nson,idtv
        real   :: tframe
c=================================================
c       end of input file
c=================================================
        real :: time
c******* Grid parameters**************************
        real :: dx2,dx3,dx1
        real :: dx2q,dx3q,dx1q
c        
        real, dimension(1:m1) :: tc,tm
        real, dimension(1:m2) :: rc,rm
        real, dimension(1:m3) :: zz,zm,g3rc,g3rm
c====================================================
c******* QUANTITIES FOR DERIVATIVES******************
        real, dimension(1:m3) :: udx3c,udx3m
c==========================================================
c******* Grid indices**************************************
        integer, dimension(1:m3) :: kmc,kpc,kmv,kpv,kup,kum
c===========================================================
c******* Metric coefficients *******************************
        real, dimension(1:m2) :: ap3j,ac3j,am3j

        real, dimension(1:m3) :: ap3ck,ac3ck,am3ck
        real, dimension(1:m3) :: ap3sk,ac3sk,am3sk
        real, dimension(1:m3) :: ap3ssk,ac3ssk,am3ssk   
c============================================================
c******* Variables for FFTW and Poisson solver****************
        real, dimension(13) :: ifx1
        real, dimension(3*m2/2+1) :: trigx1
        real, dimension(1:m2) :: ak2,ap
        real, dimension(1:m1) :: ak1,ao
        real, dimension(1:m3) :: amphk,acphk,apphk
        
c===========================================================
c******* Other variables ***********************************
        integer  :: n2m, n3m, n1m
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
        real, dimension(1:m2,1:m1) :: denbs,denbn
              
      end module param
      
c************* End of param module******************************
c===============================================================
c******* 2D arrays, dynamically allocated by each process*******
      module local_arrays
      use param
        implicit none
        real,allocatable,dimension(:,:,:) :: q1,q2,q3
        real,allocatable,dimension(:,:,:) :: pr,dens,rhs
        real,allocatable,dimension(:,:,:) :: ru1,ru2,ru3,ruro
        real,allocatable,dimension(:,:,:) :: dph,qcap,dq,hro
      end module local_arrays

c===============================================================
      module stat_arrays
       implicit none
       real,allocatable, dimension(:) :: q1_me,q1_rms 
       real,allocatable, dimension(:) :: q2_me,q3_me,q2_rms,q3_rms 
       real,allocatable, dimension(:) :: dens_me,dens_rms 
       real, allocatable,dimension(:) :: disste,dissth,densq3_me
       integer :: timeint_cdsp
      end module stat_arrays
c=====================================================       
      module stat2_param
       implicit none
       integer :: idprobe1,idprobe2
       integer :: idslab1,idslab2,idslab3,idslab4,idslab5
       integer :: jprobe1,jprobe2,jprobe3
       integer :: kslab1,kslab2,kslab3,kslab4,kslab5,kprobe1,kprobe2
       real :: probez1,probez2,prober1,prober2,prober3,slabdumpz
      end module stat2_param
c=====================================================       
      module mpih
        implicit none
        include 'mpif.h'
        integer :: ierr
        integer, parameter :: master=0
        integer :: MDP = MPI_DOUBLE_PRECISION
      end module mpih

      
      module mpi_param
        implicit none
        integer :: istart,iend,jstart,jend, kstart,kend
        integer :: jstartp,jendp
        integer :: dj,dk,mydata,mydatam
        integer :: djp
        integer, allocatable, dimension(:) :: offsetj,offsetk
        integer, allocatable, dimension(:) :: offsetjp
        integer, allocatable, dimension(:) :: countj,countk
        integer, allocatable, dimension(:) :: countjp
        integer, allocatable, dimension(:) :: countf
        integer(8), allocatable, dimension(:) :: offsetf 
      end module mpi_param
c====================================================
      module fftw_params
        use param, only: m2m,m2mh,m1m
        integer*8 fwd_plan_y, fwd_plan_z
        integer*8 bck_plan_y, bck_plan_z
        real     yr(m2m)
        complex  ya(m2mh)
        complex  zr(m1m),za(m1m)
        integer FFTW_PATIENT, FFTW_FORWARD, FFTW_BACK
        parameter (FFTW_PATIENT=32)   
        parameter (FFTW_FORWARD=-1)   
        parameter (FFTW_BACK=1)   
      end module fftw_params
