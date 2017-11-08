Matlab interfaces are now available to a growing number of GALAHAD packages. 

                            ---------
                            For LINUX
                            ---------

                       ....................
                       MATLAB before R2011a
                       ....................

To use the Matlab interfaces, GALAHAD must be installed using the
g95 compiler - other compilers may become available in the future
if Matlab's mex facility becomes more fortran-90 friendly. Either
select the "install with Matlab" option when installing, or if
a g95 version has already been installed, issue the commands

  cd $GALAHAD/src/matlab
  make -s -f $GALAHAD/makefiles/pc.lnx.g95 

(substitute pc.lnx.g95 for the appropriate string on non-Linux or 64-bit 
 machines).

N.B. the MYMATLAB environment variable must point to your system matlab 
directory; Matlab's mex executable should be found under $MYMATLAB/bin.

Once the Matlab versions have been installed, make sure that
$GALAHAD/src/matlab is on your Matlab path.

 Issue the commands

  help galahad

to find the current status of available interfaces.

For 64-bit Matlab, make sure that you are using the 64-bit g95 compiler,
and that long integers are used on the compile lines. This requires that
the -i8 and -fPIC flags be added to the BASIC variable in 
$GALAHAD/makefiles/pc64.lnx.g95, and to the FFLAGS variable in your Matlab 
mexopts.sh file (in the relevant release subdirectory of the .matlab directory 
in your home directory).

                     ........................
                     MATLAB for R2011a-R2013a
                    .........................

For MATLAB R2011a and above, GALAHAD must be installed using the
gfortran compiler. BUT ... you will need to use (and have installed)
gcc/gfortran-4.4, not the more modern version (4.5 and above) that comes 
as default with today's Linux. For Ubuntu Linux, see 

   https://help.ubuntu.com/community/MATLAB

for details on how to download packaged versions of the relevant 
outdated compilers. For other Linux distributions, you might have
to build gcc-4.4 from source code. Grumble to the Mathworks!

Once you have a working gcc/gfortran 4.4 - making sure that they are
both on your shell search path - either select the "install with Matlab" 
option when installing, or if a gfortran version has already been installed, 
issue the commands

  cd $GALAHAD/src/matlab
  make -s -f $GALAHAD/makefiles/pc.lnx.gfo 

(substitute pc.lnx.gfo for the appropriate string on non-Linux or 64-bit 
 machines).

N.B. the MATLAB environment variable must point to your system matlab directory.

Once the Matlab versions have been installed, make sure that
$GALAHAD/src/matlab is on your Matlab path.

 Issue the commands

  help galahad

to find the current status of available interfaces.


                     ...........................
                     MATLAB for R2013b and above
                    ...........................

As for MATLAB for R2011a-R2013a, but you will need gfortran/gcc-4.7 not 4.4.
Edit $GALAHAD/makefiles/pc.lnx.gfo to check that all mentions of gfortran/gcc
have the trailing -4.7.

                         ------------
                         For MAC OS X
                         ------------

Here, the supported compiler is GNU gfortran. So

  cd $GALAHAD/src/matlab
  make -s -f $GALAHAD/makefiles/pc.lnx.gfo 

will install the 32-bit version. For the 64-bit version, once again
long integers must be used, and adding -fdefault-integer-8 to the
BASIC and FFLAGS variables as described above will achieve this

                          -----------
                          For Windows
                          -----------

As we only offer support for GALAHAD in Windows using MinGW/MSYS 
GNU environent (see README.windows), and as we have not checked
whether g95/gfortran is compatible with Matlab in this case, we
are sorry but the Windows user is on her/his own. Matlab claim
to support Intel Visual Fortran as their default Windows-fortran
interface.

Nick Gould          (nick.gould@stfc.ac.uk)
Dominique Orban     (dominique.orban@polymtl.ca)
Philippe Toint      (philippe.toint@fundp.ac.be)

For GALAHAD productions
This version: 21st March 2013

