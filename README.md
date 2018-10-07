# Galahad and CUTEst test environment for optimization
Only for academic use!!!
+ https://github.com/ralna/GALAHAD/wiki
+ http://www.galahad.rl.ac.uk/

## The SIF files for test problems are downloaded from
+ svn checkout --username anonymous http://ccpforge.cse.rl.ac.uk/svn/cutest/sif/trunk ./sif
+ svn checkout --username anonymous http://ccpforge.cse.rl.ac.uk/svn/cutest/marosmeszaros/trunk ./marosmeszaros 

### Operation system: Ubuntu 14.04 with 64-bit
### Version of Matlab: R2014a
### Need: GNU gfortran 4.7, gcc 4.7

## Install the Galahad and CUTEst environment:
(1) Set the environment variable by adding the following to the file .bashrc
change the name "pirate" and Matlab path to yours

+ export ARCHDEFS="/home/pirate/archdefs"
+ export SIFDECODE="/home/pirate/sifdecode"
+ export CUTEST="/home/pirate/cutest"
+ export MASTSIF="/home/pirate/sif"
+ export PATH="${SIFDECODE}/bin:${PATH}"
+ export PATH="${CUTEST}/bin:${PATH}"
+ export MANPATH="${SIFDECODE}/man:${MANPATH}"
+ export MANPATH="${CUTEST}/man:${MANPATH}"

+ export MYARCH="pc64.lnx.gfo"
+ export MYMATLABARCH="pc64.lnx.gfo"
+ export MATLAB="/usr/local/MATLAB/R2014a"

+ export GALAHAD="/home/pirate/galahad"
+ export PATH="${GALAHAD}/bin:${PATH}"
+ export MANPATH="${GALAHAD}/man:${MANPATH}"

(2) cd $GALAHAD and run $ARCHDEFS/install_optrove

Then wait for compiling!
