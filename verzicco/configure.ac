#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([AFiD], [1.0], [john.donners@surfsara.nl])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

# AC_PROG_FC
AX_PROG_FC_MPI([true],
  [ use_mpi=yes],
  [AC_MSG_FAILURE([MPI compiler needed, but could not find MPI.],1)
  ])

AC_FC_MODULE_FLAG
# Lookup extension for modules, used for 'make clean'
AC_FC_MODULE_EXTENSION
AC_SUBST(FC_MODEXT)

# Lookup flag for passing preprocessor flags to the Fortran compiler. Used to define HAVE_CONFIG_H
AC_FC_PP_DEFINE
AC_SUBST(FC_DEFINE)

AC_PROG_AWK
AC_PROG_GREP
AC_PROG_SED

AC_LANG([Fortran])
#set the extension to f90, to compile the HDF5 test
AC_FC_SRCEXT(f90)

#Disable OpenMP as long as it is not efficient
#AC_OPENMP
#FCFLAGS="${FCFLAGS} ${OPENMP_FCFLAGS}"

#check for option to promote 4-byte REALs to 8-byte DOUBLE PRECISIONs
dnl Try these flags:
dnl   GCC >= 4.2           -fdefault-real-8 -fdefault-double-8
dnl   ifort                -r8
dnl   PGI fortran          -Mr8
dnl   xlf90                -qautodbl=dbl4
dnl   Cray ftn             -s real64
AC_MSG_CHECKING([for flag to promote REAL to DOUBLE PRECISION])
dummy_FCFLAGS=$FCFLAGS
for ax_option in -r8 -qautodbl=dbl4 "-fdefault-real-8 -fdefault-double-8" "-s real64" -Mr8; do
  FCFLAGS="$ax_option $dummy_FCFLAGS"
  AC_LINK_IFELSE(
    [AC_LANG_SOURCE([[
      program conftest
      implicit none

      print*,dsin(0.0)

      end program]])],
    [real_to_double=$ax_option])
  if test "$real_to_double" != ""; then
    AC_MSG_RESULT([$real_to_double])
    break
  fi
done

# Checks for libraries.
AX_BLAS([found_blas=yes],[AC_MSG_ERROR([Unable to find blas library.])])
AX_LAPACK([found_lapack=yes],[AC_MSG_ERROR([Unable to find lapack library.])])

AC_PATH_PROG(FFTW3_BIN,[fftw-wisdom])
# Emulating unix commands as sed: http://sed.sourceforge.net/local/docs/emulating_unix.txt
# Doubling up of [ and ] : http://lists.gnu.org/archive/html/autoconf/2011-02/msg00064.html
# original sed command: sed 's,[^/]*/[^/]*$,,'
FFTW3_ROOT=$(echo ${FFTW3_BIN}|$SED 's,[[^/]]*/[[^/]]*$,,')

AC_ARG_WITH([fftw3], AS_HELP_STRING([--with-fftw3=<root>],
  [root path to the fftw3 library]),
  [FFTW3_ROOT=$withval])
FFTW3_LDFLAGS=-L$FFTW3_ROOT/lib
AC_CHECK_LIB([fftw3], [dfftw_plan_many_dft],[FFTW3_LIBS="${FFTW3_LDFLAGS} -lfftw3"],[AC_MSG_ERROR([Unable to find fftw3 library])],[${FFTW3_LDFLAGS}])
if test "x$OPENMP_FCFLAGS" != "x"; then
  AC_CHECK_LIB([fftw3_omp], [dfftw_init_threads],[FFTW3_LIBS="${FFTW3_LDFLAGS} -lfftw3_omp -lfftw3"],[AC_MSG_ERROR([Unable to find fftw3 OpenMP-enabled library])],[${FFTW3_LIBS}])
fi
AC_SUBST(FFTW3_LIBS)

#add configure command-line argument for HDF5
AC_ARG_WITH([hdf5], AS_HELP_STRING([--with-hdf5=<root>],
  [root path to the parallel HDF5 library]),
  [HDF5_ROOT=$withval])
#look for HDF5 Fortran compiler wrapper. This is used to find installation point if HDF5_ROOT is not explicitly set.
AC_PATH_PROGS([H5FC],[h5pfc h5fc],,[${HDF5_ROOT}/bin:$PATH])
test -n "$H5FC" || AC_MSG_ERROR([Unable to find HDF5 Fortran compiler wrapper.])
# if HDF5_ROOT not set, look for it in a standard directory.
if test "x${HDF5_ROOT}" == "x"; then
  HDF5_ROOT=$(eval $H5FC -showconfig \
    | $GREP 'Installation point:' \
    | $AWK -F: '{print $[]2}' )
fi
HDF5_AM_LDFLAGS=$(eval $H5FC -showconfig \
    | $GREP 'AM_LDFLAGS:' \
    | $AWK -F: '{print $[]2}' )
HDF5_EXTRALIBS=$(eval $H5FC -showconfig \
    | $GREP 'Extra libraries:' \
    | $AWK -F: '{print $[]2}' )
HDF5_LDFLAGS="-L${HDF5_ROOT}/lib ${HDF5_AM_LDFLAGS}"

dummy_FCFLAGS=$FCFLAGS
dummy_LDFLAGS=$LDFLAGS
dummy_LIBS=$LIBS

HDF5_LIBS="-lhdf5_fortran -lhdf5 ${HDF5_EXTRALIBS}"
HDF5_FCFLAGS=${FC_MODINC}${HDF5_ROOT}/include

FCFLAGS="${HDF5_FCFLAGS}"
LDFLAGS="${HDF5_LDFLAGS}"
LIBS="${HDF5_LIBS}"

AC_MSG_CHECKING([for parallel HDF5 Fortran library])
AC_LINK_IFELSE(
  [AC_LANG_SOURCE([[
      program conftest
      use hdf5
      implicit none
      integer         :: classtype
      integer         :: comm
      integer (hid_t) :: plist_id
      integer         :: hdferr
      call h5pset_fapl_mpio_f(classtype,comm,plist_id,hdferr)
      end]])],
  [hdf5_works=yes],
  [hdf5_works=no AC_MSG_ERROR([Unable to find hdf5 Fortran library with parallel I/O])])
AC_MSG_RESULT([$hdf5_works])

FCFLAGS="${dummy_FCFLAGS} ${HDF5_FCFLAGS}"
LDFLAGS="${dummy_LDFLAGS} ${HDF5_LDFLAGS}"
LIBS="${dummy_LIBS}"
#restore FCFLAGS and add HDF5 include directory
AC_SUBST(HDF5_LIBS)

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
