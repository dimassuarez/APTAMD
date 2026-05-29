# This file is sourced by the APTAMD scripts to define environment variables 
# Location of the APTAMD scripts suite
APTAMD="/home/studgeuo/APTAMD/"
# AMBERTOOLS and AMBER packages compiled with MPI, NETCDF, OPENMP, and CUDA options 
# We are AMBER fans: https://ambermd.org
AMBERHOME="/opt/apps/AL9/amber26"
# Parallel BASH is needed to distribute tasks among available procesors. 
# https://www.gnu.org/software/parallel/sphinx.html
PARHOME="/opt/apps/AL9/parallel-bash"
# Scratch space. Many temporary directories and files are generated
SCRATCH="/scratch"
# The TOOLS directory contains several in-house fortran codes that perform auxiliary tasks
# hopefully they will be replaced by Python scripts soon
TOOLS="$APTAMD/AUXTOOLS"
# Some structural analysis of DNA fragments are performed using the X3DNA software
# https://x3dna.org/
X3DNAHOME="/opt/apps/AL9/x3dna-v2.6/"
# DSSP  https://github.com/PDB-REDO/dssp: not needed for aptamers
# DSSPHOME="/opt/dssp"

# OpenMPI (needed for parallel execution of AMBER programs)
# AMBER may also include OpenMPI 
# https://www.open-mpi.org/
# Select here the same OpenMPI used to compile AMBER programs!
MPI_HOME=/usr/lib64/openmpi
LD_LIBRARY_PATH=$MPI_HOME/lib:$LD_LIBRARY_PATH

# Conformational entropy: cencalc entropy software  
# (pulled automatically during installation)
CENCALC_PATH=$APTAMD/AUXTOOLS/cencalc_quicksort

# OCTAVE  Most of the data analysis and plotting are performed using octave scripts
# https://octave.org/
# At some point all the octave code in APTAMD will be replaced by python 
OCTAVE=$(which octave | grep -v alias)
# Some data handling is performed using datamash
DATAMASH=$(which datamash | grep -v alias)

#  Docking calcs are handled using the Autodock suite and 
#  the Autodock tools included in mgltools 
# https://autodock.scripps.edu
ADCK="/opt/apps/AL9/adck_tools/"
MGLTOOLS="/opt/apps/AL9/mgltools/"
#
# Python dependencies:
# GAMD reweight as well as parmed require PYTHON with numpy and matplotlib 
# GAMD python scripts work with either version 2 or 3)
# TSNE clustering requires python3 with numpy, matplotlib, netcdf4, scikit-learn 
# It is recommended to use the miniconda distro in AMBER with augmented packages (scikit and netcdf4) 
# Be careful with the PYTHONPATH variable
export PYTHON=$AMBERHOME/miniconda/bin
export PYTHONPATH=$AMBERHOME/lib/python3.12/site-packages/

# For MMPBSA calcs. 
# MSMS program is used in the non-polar components of MMPBSA calcs 
# https://ccsb.scripps.edu/msms/
MSMS=/opt/apps/AL9/msms/msms

# For QMMMPBSA calcs (only). 
# https://orcaforum.kofo.mpg.de/app.php/portal
# ORCA=/opt/apps/AL8/orca6
# D3H4 corrections in SCC_DFTB3 require CUBY4  http://cuby4.molecular.cz/
# CUBY4=/opt/cuby4/cuby4

if [   -z "$AMBERHOME" ]; then echo 'AMBER is not available, but needed!'; exit; fi
if [   -z "$OCTAVE"    ]; then echo 'OCTAVE is not available, but needed!'; exit; fi
if [ ! -e "$SCRATCH"   ]; then echo '$SCRATCH space does not exist, but needed!'; exit; fi

