# This file is sourced by the APTAMD scripts to define environment variables 

# Location of the APTAMD scripts suite
APTAMD="/home/dimas/SCRIPTS/SUITE_APTAMD/"
# AMBERTOOLS and AMBER packages compiled with MPI, NETCDF, OPENMP, and CUDA options 
# We are AMBER addicts:  https://ambermd.org
AMBERHOME="/opt/amber23"
# Parallel BASH is needed to distribute tasks among available procesors. 
# https://www.gnu.org/software/parallel/sphinx.html
PARHOME="/opt/parallel-bash/"
# Scratch space. Many temporary directories and files are generated
if [ -n "$PBS_ENVIRONMENT" ] ; then
  echo "scratch=$SCRATCH"
else
  SCRATCH="/scratch"
fi
# The TOOLS directory contains several in-house fortran codes that perform auxiliary tasks
TOOLS="$APTAMD/AUXTOOLS"
# Docking cacls are handled using the Autodock suite and Autodock tools included in MGLTOOLS
# https://autodock.scripps.edu/download-autodock4/
# https://ccsb.scripps.edu/mgltools/
ADCK="/opt/autodock/"
MGLTOOLS="$APTAMD/mgltools_x86_64Linux2_1.5.7/"
# Some structural analysis of DNA fragments are performed using the X3DNA software
# x3dna is now a commercial program. 
# https://x3dna.org/
X3DNAHOME="$APTAMD/x3dna-v2.4/"
# DSSP  https://github.com/PDB-REDO/dssp
DSSPHOME="/opt/dssp"
# OpenMPI (needed for parallel execution of AMBER programs)
# AMBER may also include OpenMPI 
# For example https://www.open-mpi.org/
MPI_HOME=/opt/openmpi/4.1.1-gcc85/
LD_LIBRARY_PATH=$MPI_HOME/lib:$LD_LIBRARY_PATH
# Most of the data analysis and plotting is performed using octave scripts
# https://octave.org/
# At some point all the Octave code in APTAMD will be replaced by python 
OCTAVE=$(which octave | grep -v alias)
# Some data handling is performed using datamash
DATAMASH=$(which datamash | grep -v alias)
# GAMD reweight as wellas parmed require PYTHON with numpy and matplotlib 
# GAMD python scripts work with either version 2 or 3)
# TSNE clustering requires python3 with numpy, matplotlib, netcdf4, scikit-learn 
# Better use AMBER Python with augmented packages (scikit and netcdf4) 
export PYTHON=$AMBERHOME/miniconda/bin
export PYTHONPATH=$AMBERHOME/lib/python3.10/site-packages/

# Conformational ENTROPY calculations https://github.com/dimassuarez/cencalc_quicksort
CENCALC_PATH="/home/dimas/SCRIPTS/CCMLA_QSORT/"
DISLIN="/opt/dislin/" # https://www.dislin.de/ used by cencalc for plotting 

# For QMMMPBSA and MMPBSA calcs:
# ORCA (for QM/MM calcs using do_mmpbsa) https://orcaforum.kofo.mpg.de/app.php/portal
ORCA=/opt/orca5
# Terachem (if available may be used for QM/MM calcs) http://www.petachem.com/products.html
TeraChem=/opt/TeraChem
# MSMS program is used in the MMPBSA calcs https://ccsb.scripps.edu/msms/(
MSMS=/home/dimas/SCRIPTS/APTAMD/MSMS/msms
# D3H4 corrections in SCC_DFTB3 require CUBY4  http://cuby4.molecular.cz/
CUBY4=/opt/cuby4/cuby4

if [   -z "$AMBERHOME" ]; then echo 'AMBER is not available, but needed!'; exit; fi
if [   -z "$OCTAVE"    ]; then echo 'OCTAVE is not available, but needed!'; exit; fi
if [ ! -e "$SCRATCH"   ]; then echo '$SCRATCH space does not exist, but needed!'; exit; fi

