The FORTRAN codes in this directory implement some data transformations ,
geometric analysis, empirical calculations , etc.   that are requested by
the driver scripts in the MMPBSA, ENTROPY and DOCKING directories. Of course,
these are indeed outdated codes and will be replaced by other tools in 
(hopefully) future revisions of the suite. Dynamic memory allocation
is not used in most of these codes, but the default sizes of the parameter
variables should be large enough for most practical purposes.

Execute the comp.sh to compile these auxiliary programs with the GNU gfortran compiler

./comp.sh all
