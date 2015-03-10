#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([AFiD], [1.0], [john.donners@surfsara.nl])
AC_CONFIG_HEADERS
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE

# AC_PROG_FC
AX_PROG_FC_MPI([true],
  [ use_mpi=yes],
  [AC_MSG_FAILURE([MPI compiler needed, but could not find MPI.],1)
  ])

AC_FC_MODULE_FLAG

AC_PROG_AWK
AC_PROG_GREP

AC_LANG([Fortran])
AC_OPENMP
FCFLAGS="${FCFLAGS} ${OPENMP_FCFLAGS}"


# Checks for libraries.
AX_BLAS([found_blas=yes],[AC_MSG_ERROR([Unable to find blas library.])])
AX_LAPACK([found_lapack=yes],[AC_MSG_ERROR([Unable to find lapack library.])])

# check for fftw3 library with C-compiler. Fortran would look for symbols with underscores
AC_LANG([C])

AC_ARG_WITH([fftw3], AS_HELP_STRING([--with-fftw3=<root>],
  [root path to the fftw3 library]),
  [FFTW3_ROOT=$withval FFTW3_LDFLAGS=-L$FFTW3_ROOT/lib])
AC_CHECK_LIB([fftw3], [fftw_plan_guru_dft],[FFTW3_LIBS="${FFTW3_LDFLAGS} -lfftw3"],[AC_MSG_ERROR([Unable to find fftw3 library])],[${FFTW3_LDFLAGS}])
if test "x$OPENMP_FCFLAGS" != "x"; then
  AC_CHECK_LIB([fftw3_omp], [fftw_init_threads],[FFTW3_LIBS="${FFTW3_LDFLAGS} -lfftw3_omp -lfftw3"],[AC_MSG_ERROR([Unable to find fftw3 OpenMP-enabled library])],[${FFTW3_LIBS}])
fi
AC_SUBST(FFTW3_LIBS)

AC_ARG_WITH([hdf5], AS_HELP_STRING([--with-hdf5=<root>],
  [root path to the parallel HDF5 library]),
  [HDF5_ROOT=$withval HDF5_LDFLAGS=-L${HDF5_ROOT}/lib])
AC_PATH_PROGS([H5FC],[h5pfc h5fc],[${HDF5_ROOT}/bin:$PATH])
# if HDF5_ROOT not set, look for it in a standard directory.
if test "x${HDF5_ROOT}" == "x"; then
  HDF5_ROOT=$(eval $H5FC -showconfig \
    | $GREP 'Installation point:' \
    | $AWK -F: '{print $[]2}' )
  HDF5_LDFLAGS=-L${HDF5_ROOT}/lib
fi

FCFLAGS="${FCFLAGS} ${FC_MODINC}${HDF5_ROOT}/include"
AC_CHECK_LIB([hdf5_fortran], [h5pset_fapl_mpio_c_],[HDF5_LIBS="${HDF5_LDFLAGS} -lhdf5_fortran -lhdf5"],[AC_MSG_ERROR([Unable to find hdf5 Fortran library with parallel I/O])],[${HDF5_LDFLAGS} -lhdf5])
AC_SUBST(HDF5_LIBS)

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
