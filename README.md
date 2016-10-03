AFiD
====
AFiD is a highly parallel application for Rayleigh-Benard and Taylor-Couette flows. 
See [van der Poel et al (2015)](http://dx.doi.org/10.1016/j.compfluid.2015.04.007)
for more details. 

It is developed by Twente University, SURFsara and University of Rome "Tor Vergata". 

See the file COPYING for copying permission of the 2DECOMP&FFT library.
This application is free and unencumbered software released into the public domain. 

Installation
------------
The AFiD model has the following prerequisites:

 * MPI
 * BLAS
 * LAPACK
 * FFTW3
 * HDF5 with parallel I/O 

It's recommended to download a release tarball of AFiD, which can be found ​[here](https://github.com/jdonners/afid/releases). To install AFiD, please 
use the 'configure' script. Note that you'll need to set optimization and debugging options yourself. 
The easiest way to configure and build AFiD is

```
./configure
make
make install prefix=/path/to/install/afid
```

It tries to find and configure all prerequisites automatically, although it doesn't always succeed. 
By default it uses the -O2 optimization flag (if available).
The most important configuration options are:

```
./configure MPIFC=mpif90.gfortran              # set MPIFC to your MPI compiler wrapper for Fortran
./configure --with-blas=/path/to/blas.lib      # library with blas routines
./configure --with-lapack=/path/to/lapack.lib  # library with lapack routines 
./configure FCFLAGS=-O3                        # very high optimization
./configure FCFLAGS="-g -O0"                   # debug info, no optimization
```

The configure script locates the fftw-wisdom utility to find the root path of the FFTW3 library and it uses the h5pfc compiler wrapper 
to configure the HDF5 library. You can override these using:

```
./configure --with-fftw3=<root path to fftw3 installation>
./configure --with-hdf5=<root path to hdf5 installation> 
```

It is recommended to use the vendor-optimized libraries for BLAS and (possibly) LAPACK (e.g. MKL, ESSL or LibSci). 
Note that the FFTW3 library cannot be replaced with the MKL library, since it doesn't support the calls that are used in AFiD.

Should you want to build from the repository (when you download the source file from directly using "Download ZIP") you first need to create the configure script. For this you need a recent versions of the GNU autotools and the configure script is created using the command 

```
autoreconf
```

Usage
-----
See the [manual](MANUAL.md) for more a description of the input parameters of the code. An example of the usage of the code can be found ​[here](https://github.com/jdonners/afid/tree/master/Example). 
