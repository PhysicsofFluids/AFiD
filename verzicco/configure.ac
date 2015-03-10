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

AC_PATH_PROG([H5CC],[h5cc])
HDF5_tmp_inst=$(eval $H5CC -showconfig \
   | $GREP 'Installation point:' \
   | $AWK -F: '{print $[]2}' )
FCFLAGS="${FCFLAGS} ${FC_MODINC}${HDF5_tmp_inst}/include"

AC_OPENMP
AC_LANG([Fortran])
AC_OPENMP

# Checks for libraries.
AX_BLAS([found_blas=yes],[AC_MSG_ERROR([Unable to find blas library.])])
AX_LAPACK([found_lapack=yes],[AC_MSG_ERROR([Unable to find lapack library.])])

# check for fftw3 library with C-compiler. Fortran would look for symbols with underscores
AC_LANG([C])

AC_CHECK_LIB([fftw3], [fftw_plan_guru_dft],[FFTW3_LIBS="-lfftw3"],[AC_MSG_ERROR([Unable to find fftw3 library])])
if test "x$OPENMP_CFLAGS" != "x"; then
  AC_CHECK_LIB([fftw3_omp], [fftw_init_threads],[FFTW3_LIBS="-lfftw3_omp -lfftw3"],[AC_MSG_ERROR([Unable to find fftw3 OpenMP-enabled library])],[-lfftw3])
fi
AC_SUBST(FFTW3_LIBS)

AC_CHECK_LIB([hdf5_fortran], [h5pset_fapl_mpio_c_],[HDF5_LIBS="-lhdf5_fortran -lhdf5"],[AC_MSG_ERROR([Unable to find hdf5 Fortran library with parallel I/O])],[-lhdf5])
AC_SUBST(HDF5_LIBS)

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
