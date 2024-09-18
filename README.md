# **APTAMD**

### **Description** 
The application of a computational protocol aimed to refine the 3D structure of aptamer molecules is not straightforward because multiple programs and a large variety of datasets are used. APTAMD is a set of specific scripts that are designed to automate and streamline the most important stages of an MD-based aptamer building protocol that relies on the [AMBER](https://ambermd.org/) suite of programs: (1) building of an initial model starting from a preliminary 3D structure, (2) Gaussian-accelerated MD simulation, (3) structural analysis and construction of a free energy map, (4) conventional MD simulation, and (5) clustering and energy scoring of the MD trajectory of the aptamer. 

<img src="./Images/APTAMD_Protocol.png" width="800" height="450" style="display: block; margin: 0 auto">
 
The computational protocol implemented in APTAMD is described in detail in the following references:

* A. Díaz-Fernández, R. Miranda-Castro, N. Díaz, D. Suárez, N. de-los-Santos-Álvarez and M.J. Lobo-Castañón. (2020). Aptamers targeting protein-specific glycosylation in tumor biomarkers: general selection, characterization and structural modeling. Chemical Science 11, 9402-9413. [DOI](https://doi.org/10.1039/D0SC00209G )

* A. Díaz-Fernández, C. S. Ciudad, N. Díaz, D. Suárez, N. de-los-Santos-Álvarez and M.J. Lobo-Castañón. (2024). Refinement and Truncation of DNA Aptamers based on Molecular Dynamics Simulations: Computational Protocol and Experimental Validation. *Submitted*
 
All questions regarding the usage of APTAMD or bug reports should be addressed to Dimas Suárez (dimas@uniovi.es).

### **Installation**

The suite is distributed as collection of scripts for Linux systems. Most of the scripts are written in [BASH](https://www.gnu.org/software/bash/) altough some numerical and plotting tasks are performed using Octave or Python scripts. There are also a few Fortran codes that carry out auxiliary tasks. In a few months, we will replace all the Fortran and Octave codes by Python scripts, but the BASH backbone will be maintained.    

To install APTAMD, download the suite with the commands:

`git clone https://github.com/dimassuarez/APTAMD.git`
   
Define the `APTAMD` environmental variable pointing to the directory containing the ATPAMD files:

`export APTAMD=/mydir/APTAMD` 

Compile the auxiliary Fortran codes: 

`cd $APTAMD/AUXTOOLS; ./comp.sh all`

Edit the `$ATPAMD/ENV/aptamd_env.sh` file and adjust the BASH variables pointing to the software tools used by APTAMD (Of course AMBER and other software have to be installed on your system). 

Add `$ATPAMD/bin` to your `$PATH` environment variable:

`export PATH=$APTAMD/bin:$PATH`

Entering at the command line the name(s) of the APTAMD scripts (e.g., `do_aptamer_edition`) should print out a brief description of their functions. 

**CONTENTS**

The APTAMD collection is organized in the following subdirectories:

* EDITION           : Scripts to prepare coordinate and parameter files 
* RUNMD             : Scripts and files to perform GaMD and MD simulations 
* STRUCT            : Scripts to extract MD frames and/or carry out structural analysis and clustering  
* RWGAMD            : Scripts and files to carry out the energy reweighing of the GaMD trajectories
* MMPBSA            : Scripts to perform MMPBSA-like calculations
* AUXTOOLS          : Scripts and auxiliary Fortran codes.
* ENV               : aptamd_env.sh is located here
* bin               : Links to the main BASH scripts
* DOC               : ASCII text files providing help info and details of the computational protocol

**DEPENDENCIES**

* Compiler and Script Interpreters: Octave (version >= 5) , Python3, GNU Parallel Bash, GNU GFortran
* AMBER Suite
* X3DNA-DSSR
* MSMS 

Note: The MMPBSA scripts can handle QM/MM calculations using various QM codes (see comments in [aptamd_env.sh](ENV/aptamd_env.sh) ). 

### **Hardware and OS requirements**
The recommended starting configuration for a workstation or a computer server dedicated to perform molecular simulations of systems with the typical size of aptamers, would include a dual socket motherboard equipped with 2 multicore CPUs, 4 GB of RAM per CPU core and a storage capacity greater than 4 TB constituted by SSD/NVMe drives. The same computer should be equipped with at least two state-of-the-art NVIDIA GPU cards. Assuming that a Linux operating system (e.g. [Almalinux](https://almalinux.org/)) and the AMBER package are installed, such a workstation/server machine would be able to simultaneously carry out production runs on the GPUs with minimum CPU usage and other tasks for the preparation and analysis of the GaMD/MD trajectories using the remaining RAM and CPU cores. Further information pertaining to the recommended hardware can be found on the AMBER website [GPU info](https://ambermd.org/GPUHardware.php). All major HPC hardware vendors offer workstations and clusters suitable for performing MD simulations.

---



## **User Guide**

### Initial model: `do_aptamer_edition` 
For a given aptamer sequence, its secondary (2D) structure can be predicted using the [mfold](http://www.unafold.org/DNA_form.php) algorithm. From the 2D mfold structure, initial 3D coordinates in PDB format can be obtained using the [RNA Composer webserver](https://rnacomposer.cs.put.poznan.pl/). The preliminary 3D models are processed automatically using `do_aptamer_edition.sh` that executes and monitors other scripts and programs of the [AMBER suite](https://ambermd.org/). `do_aptamer_edition.sh` transforms the atom and residue names of RNA to those of DNA, removes all 2’-hydroxyl groups on the ribose sugars, adds the missing methylene in thymines and all the H atoms, and assigns the required Molecular Mechanics (MM) parameters from the [parmbsc1 force field](https://mmb.irbbarcelona.org/www/ParmBSC1),  relaxes the internal geometry of the nucleobases and adds an octahedral box  water molecules including and Na+/Cl- counterions. 



---

<div style="background-color: rgb(220,220,220);">

**EXAMPLE**

In this example, we build a *denovo* model for the unbound form of a truncated sequence corresponding to the [anti-MUC1 DNA aptamer](https://doi.org/10.1111/j.1742-4658.2011.08440.x).  In particular, we adopt the following 26-mer sequence: `3'-GCAGTTGATCCTTTGGATACCCTGGT-5'`. The [2L5K](https://www.rcsb.org/structure/2l5k) structure in the PDB contains several NMR models of a closely-related sequence. Both input and output files generated with the APTAMD tools for this example can be accssed at the Mendeley data repository doi: 10.17632/hg98523cbs.1

The first step is to enter the 26-mer sequence in the [mfold server](http://www.unafold.org/mfold/applications/dna-folding-form.php) selecting T=25 oC and 0.150 M NaCl ionic strength. Leave other settings at their defaults.  The following 2D structure should be obtained.

<img src="./Images/2l5k_mfold.png" width="200" height="100" title="Mfold 2d model" style="display: block; margin: 0 auto" >

From the mfold output, we copy the `.ct` file and paste it into the box of the Format converter (CT to dot-bracket) of the [RnaComposer server](https://rnacomposer.cs.put.poznan.pl/tools) to obtain the following secondary (2D) structure in dot-bracket notation:    

<section style="font-family:'Courier New'">
GCAGTTGATCCTTTGGATACCCTGGT
    
.......((((...))))........
</section>

The 2D structure is transformed into a 3D RNA model. Copy and paste it in the last two lines of the interactive mode box. After a few seconds, a first 3D structure is generated. Save the initial coordinates of the model in a PDB file. 

<img src="./Images/2l5k_rnacomposer.png" width="300" height="150" title="RNAComposer 3D model" style="display: block; margin: 0 auto"  >


Select or make a working directory that will contain the main directory of the simulation to be done. Use then the following commands:

`cd` *Working Directory*

`mkdir 2L5K_model_GAMD; cd 2L5K_model_GAMD; mkdir 1.EDITION; cd 1.EDITION`
`cp $APTAMD/EXAMPLE_INPUT_FILES/2L5K_model_initial.pdb .`

where `2L5K_model_initial.pdb` is the PDB file with the initial coordinates generated by RNA composer using the 2D mfold structure (we're not using the PDB data here!). 

Assuming that `$APTAMD/bin` is in your PATH, the MM parameter files are obtained with:

`do_aptamer_edition  2L5K_model_initial.pdb`

what generates multiple files:

<section style="font-family:'Courier New'">
2L5K_model_initial.pdb  2L5K_model_solute.pdb   edit_leap_solute.src     sander_relax_solute. 
2L5K_model.pdb          2L5K_model_solute.top  edit_leap.log         edit_leap.src
2L5K_model.crd  2L5K_model_solute.crd   2L5K_model.top , ....
</section>


For example, `2L5K_model.top` and `2L5K_model.crd` are the topology (parameter) and intial coordinate files of the solvated aptamer ready to be used in the subsequent GaMD simulation.

</div>



### Gaussian Accelerated MD simulation : `do_runmd` 
The conformational space of the solvated aptamer is explored by means of a Gaussian Accelerated Molecular Dynamics (GaMD) simulation or by conventional MD simulations. GaMD uses harmonic boost potentials to smooth out the potential energy surface, accelerating thus transitions between low-energy configurations. `do_runmd.sh` controls various preparatory stages (solvent relaxation, thermalization, pressurization and GaMD equilibration) and finally launches the production phase of the GaMD/MD simulations. The settings of the GaMD and MD calculations driven by do_runmd.sh are normally employed in many biomolecular simulations. See [here](SUITE_APTAMD/DOC/MD_settings.txt) for further details.  

---
<div style="background-color: rgb(220,220,220);">
    
**EXAMPLE**

To launch the GaMD simulation of the 2L5K model

`cd` *Working directory*    

`cp $APTAMD/EXAMPLE_INPUT_FILES/input_rungamd.src  2L5K_model_GAMD/`

The basic options for this GaMD job are described in the [input_rungamd.src file](EXAMPLE_INPUT_FILES/input_rungamd.src). Then we run the job: 

`do_runmd  2L5K_model_GAMD/input_rungamd.src`

*WARNING* 1-2 days of GPU usage and multicore CPU are required to run the GAMD simulation.  

`do_runmd` makes new subdirectories into the project directory (i.e. 2L5K_model_GAMD) and writes many input and output files. Most importantly, the `5.PRODUCTION` folder contains the trajectory files named as md_001.* , md_002.* etc. For example, 

<section style="font-family:'Courier New'">
gamd_eq.mdcrd    ---> GaMD equilibration coordinates <br />
gamd_eq.out      ---> GaMD equilibiration output<br />
gamd.log<br />
gamd-restart.dat<br />
gamd.rst<br />
job_GAMD.sh<br />
md_000.rst      ---> Starting point of the production phase = gamd.rst <br />
md_001.gamd_log ---> Log file with GaMD boost potential data<br />
md_001.out      ---> PMEMD output files <br />
md_001.rst      ---> PMEMD restart file<br />
md_001_solute.mdcrd     ---> Coordinates of the solute atoms along the simulation  <br />
md_001_solutewat.mdcrd  ---> Coordinates of solute, counterions and a shell of waters <br />
md_002.gamd_log<br />
md_002.out<br />
md_002.rst<br />
md_002_solute.mdcrd<br />
md_002_solutewat.mdcrd<br />
...
</section>



The folders 2.RELAX_SOLVENT, 3.THERMALIZATION and 4.PRESSURIZATION contain the intermediate files produced during the preparatory stages of the GaMD simulation.  

*Hint:* If necessary the GaMD simulation can be resumed by changing & executing the `job_GaMD.sh` script in 5.PRODUCTION.

*Warning:* Due to the stochastic control of temperature (Langevin's method) exerted during the GaMD/MD simulations, two MD simulations started from identical coordinates will be different. Of course average properties should be comparable.   
</div>


### Analysis of GaMD trajectories: `do_struct`  and `do_reweight_gamd`
The GaMD conformations can be characterized by two structural indexes: 

* the root-mean-squared-deviation (RMSD) of the heavy atoms (P, C, N, O) with respect to the initial structure.
* the interaction network fidelity index (INF), which is built from the sets of characteristic intramolecular interactions in a reference structure ($S_r$) and in a given snapshot ($S_m$) (See [Parisien et. al.](https://doi.org/10.1261%2Frna.1700409)). The values of $S_r$ and $S_m$, are determined by the [DSSR](https://x3dna.org/) software, which identifies both base pair interactions and non-pair interactions (i.e., base stacking or base contacts).

A single script `do_struct.sh`  processes all the GaMD trajectory files, evaluates these and other structural descriptors using  the cpptraj and DSSR programs, calculates the statistical average and standard deviations of the RMSD/INF data and shows graphically their evolution along the GaMD trajectory.

---
<div style="background-color: rgb(220,220,220);">
    
**EXAMPLE**
To obtain the RMSD plot , the RGYR plot (radius of gyration) and the INF plots of the aptamer during the GaMD simulation:

`cd` *Working directory*    

`cp $APTAMD/EXAMPLE_INPUT_FILES/input_struct.src  2L5K_model_GAMD/`

`do_struct  2L5K_model_GAMD/input_struct.src`

**WARNING** The DSSR program, which is needed to calculate the DNA structural descriptors, demands much CPU time. 

Some  options for `do_struct`are specified in [input_struct.src](EXAMPLE_INPUT_FILES/input_struct.src) and the output is written in a new folder 6.ANALYSIS/STRUCT. The evolution of the RMSD/RGYR/INF descriptors is plotted along the GaMD trajectory (see the .png files). 

<img src="./Images/struct.png" width="800" height="450"  style="display: block; margin: 0 auto">

</div>

---

The GaMD simulation is energetically reweighed in terms of the RMSD/INF coordinates to produce a 2D free energy map of the conformational space. The `do_reweight_gamd.sh`script  assembles all the necessary data files, selects a proper number of bins along the RMSD/INF coordinates and plots the free energy ($G$) and the logarithm of the unweighted population ($\log_{10}(P^*)$) distributions over the 2D bins.  `do_reweight_gamd.sh` also selects a set of representative structures from the GaMD trajectory files and saves them in PDB format . 

---

**EXAMPLE**
<div style="background-color: rgb(220,220,220);">
Processing of the GaMD trajectory can be done in just one step:

`cd` *Working directory*    

`cp $APTAMD/EXAMPLE_INPUT_FILES/input_rwgamd.src  2L5K_model_GAMD/`

`do_reweight_gamd  2L5K_model_GAMD/input_rwgamd.src`

The output from `do_reweight_gamd` is saved into 6.ANALYSIS/RW_GAMD_ENE. It is normal to try different settings for reweighing (e.g., bin size or descriptors) so that `do_reweight_gamd` creates specific directories (see also an example of [input_rwgamd.src](EXAMPLE_INPUT_FILES/input_rwgamd.src)). Again many output files are produced, but 2D_RW.dat deserves particular attention because it identifies the most-likely GaMD structures using a syntax as:  
<section style="font-family:'Courier New'">
 Picking up snapshots bin= 10 for  E= 1.525000 POP= 0.030  X=  0.52  Y=  7.54 <br />
   Snapshot 1   isnap= 33922   ID=33922  RC1= 0.522976 RC2= 7.681600 Boost= 2.736762<br /> 
   bin= 10  isnap_in_crd= 13922  CRD_FILE= 2 <br />
   ... <br />
</section>

<img src="./Images/rwgamd.png" width="800" height="450" style="display: block; margin: 0 auto">

In the example, snapshot # 33922 assigned to bin #10 (located on the free energy minimum; bins with marginal populations on the biased distribution are ignored) corresponds to snapshot #13922 in the MD trajectory segment #2. It has a low value of GaMD boost potential and could be selected for subsequent cMD jobs. The coordinates are saved in a separate PDB file named as `snap_10_13922.pdb`. Note that 

*Selected GAMD model*

<img src="./Images/2l5k_gamd_model.png" width="300" height="100" style="display: block; margin: 0 auto">

</div>


###  Equilibrium properties of aptamer models from conventional MD: `do_cluster` and `do_mmpbsa` 

Conventional MD simulation (cMD) provide (valuable) equilibrium conformational sampling of aptamer models. After the most likely conformer(s) generated by  GaMD  is processed by `do_aptamer_edition.sh`. use again `do_runmd.sh`to drive a cMD simulation that usually extends up to several  µs   to ensure that the aptamer molecule relaxes and explores its equilibrium conformations in aqueous solution. 

---

<div style="background-color: rgb(220,220,220);">

**EXAMPLE**
To start a cMD from a proper GaMD structure, it is necessary to make a new "main directory" and repeat the preparatory steps like:

cd` *Working Directory*

`mkdir 2L5K_model_MD; cd 2L5K_model_MD; mkdir 1.EDITION; cd 1.EDITION`
`cp $APTAMD/EXAMPLE_INPUT_FILES/selected_file.pdb 2L5K_model_initial.pdb `

where `selected file.pdb` is the selected GaMD structure (only solute atoms).  A new edition step is then executed to build the topology files and the corresponding solvent box:

`do_aptamer_edition  2L5K_model_initial.pdb`

The conventional MD can be launched with:

`cd` *Working directory*    

`cp $APTAMD/EXAMPLE_INPUT_FILES/input_runmd.src  2L5K_model_MD/`

`do_runmd  2L5K_model_MD/input_md.src`

The variable TYPE in [input_md.src](EXAMPLE_INPUT_FILES/input_md.src) is now declared as TYPE=MD.  

*WARNING* 3-4 days of GPU usage and multicore CPU are required to run the MD simulation.  

</div>
---

Using `do_struct.sh` and `do_mmpbsa.sh` the evolution of structural and/or energetic properties are monitored along the cMD simulation. Eventually, these descriptors may exhibit pronounced drifts at the beginning of the cMD  and, therefore,  equilibrium properties and clustering analysis would be better evaluated using only the  “fully-relaxed” part of the trajectory.

---
<div style="background-color: rgb(220,220,220);">
**EXAMPLE**
Once that the cMD is completed, you can use (and adapt) the following script input files for analysis :  

`cd` *Working directory*    

`cp $APTAMD/EXAMPLE_INPUT_FILES/input_struct_md.src  2L5K_model_MD/`

`cp $APTAMD/EXAMPLE_INPUT_FILES/input_snap.src       2L5K_model_MD/`

`cp $APTAMD/EXAMPLE_INPUT_FILES/input_mmpbsa.src     2L5K_model_MD/`

The RMSD and RGYR values are calculated, but we skip now the INF ones (see for example [input_struct_md.src](EXAMPLE_INPUT_FILES/input_struct_md.src) ). 

`do_struct 2L5K_model_MD/input_struct_md.src`

Note that the resulting plot suggests a structural transition during the cMD.

<img src="./Images/struct_md.png" width="400" height="225" style="display: block; margin: 0 auto">

Before carrying out the MM-PBSA calculations, MD snapshots are to be extracted from the trajectory files using the `do_snapshots`script and saved in PDB format. 

`do_snapshots 2L5K_model_MD/input_snap.src`vi 

The PDB files for the snapshots are saved in `6.ANALYSIS/SNAPSHOTS`. In the same directory, other files contain data characterizing the Na+···APT contacts that can be of particular interest for aptamer simulations.  

Next `do_mmpbsa.sh` takes coordinates from the snapshot PDBs to calculate MM-PBSA like energies (many options can be selected, see [do_mmpbsa.txt](DOC/do_mmpbsa.txt). 

`do_mmpbsa 2L5K_model_MD/input_mmpbsa.src`

SANDER/PBSA output files (packed into a `OUTPUT.tar` file) and many other files with parsed data are saved in the 6.ANALYSIS/MMPBSA directory. Since different inner dielectric constants can be used (default is PDIE=1), the MMPBSA files are actually written in subdirectories named as PDIE_1, PDIE_4, etc.  The '.dat' files contain the data calculated for each MD snapshot while those files with the `.med` extension contain the mean values and statistical uncertainties (see the example below).  Statistics can be redone for a particular segment of trajectory using the `$APTAMD/MMPBSA/stat_plot.sh` script. See the following G_MMPBSA.med example:

<section style="font-family:'Courier New'">
 PREFIX G_MMPBSA(CMPLX) USING ALL DATA NDAT= 400<br />
 Mean       -4614.189  <br />
 Max        -4541.095  <br />
 Min        -4719.315  <br />
 SE             1.409 <br />
 BE             0.000  <br />
 LBE            2.403  <br />
#SE: standard error <br />
#BE: block error estimate (non-linear fitting) <br />
#LBE: limiting block error estimate <br />
</section>

The corresponding `.png` files display the time evolution of the various energy components.  For example,   

<img src="./Images/G_MMPBSA.png" width="400" height="225" style="display: block; margin: 0 auto">
</div>
--- 

#### About MM-PBSA

Many variants of the [MM-PBSA](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.8b00805) approach are routinely applied in a broad range of biomolecular modeling applications and there exist several tools for streamlining this type of calculations on MD ensembles. In particular, the scoring of the aptamer models is readily available by computing the conventional MM-PBSA energy of the solute atoms as:

$$ G = E_{MM} + \Delta G_{solv}^{PB} +  \Delta G_{solv} ^{np} $$

where  $E_{MM}$ is the molecular mechanics energy including the 3RT contribution due to six translational and rotational degrees of freedom, $ \Delta G_{solv}^{PB} $  is the electrostatic solvation energy obtained from Poisson-Boltzmann calculations,and $\Delta G_{solv} ^{np}$  is the non-polar part of solvation energy due to cavity formation and dispersion interactions between the solute and the solvent molecules.

---

#### Clustering

Clustering analysis of the cMD trajectory can yield critical information to optimize the sequence of the aptamer. This task can be performed using the average-linkage clustering algorithm as implemented in the [cpptraj](https://amberhub.chpc.utah.edu/cpptraj/) program. Our script  `do_cluster.sh` organizes all the input/output data files and assigns optimal values to the clustering options. The distance metric between MD frames is calculated via best-fit coordinate RMSD using the coordinates of the heavy atoms (C, N, O, P). To select a proper RMSD threshold, it is useful to test several values around 1-5 Å and choose the threshold that gives a small number of populated clusters.

---

**EXAMPLE**
<div style="background-color: rgb(220,220,220);">

`cd` *Working directory*    

`cp $APTAMD/EXAMPLE_INPUT_FILES/input_cluster.src  2L5K_model_MD/`

`do_cluster  2L5K_model_MD/input_cluster.src`

The variable EPS in [input_cluster.src](EXAMPLE_INPUT_FILES/input_cluster.src) is declared as an array (EPS="3.0 4.0 5.0") so that three different thresholds are tested using only the last half of the cMD trajectory. The aptamer model has a significant flexibility and the EPS=5.0 threshold is appropriate to yield a few populated clusters (see the summary.dat file). The cluster representatives in PDB format can be readily visualized using molecular graphics software such as Pymol or Chimera.  For example, the superposition of the top 3 clusters that account for 82% of the MD frames shows that the 3'-terminal residues adopt different conformations. When superposing the APTAMD cluster representativees, the RNAComposer model and the NMR  experimental model, we also see that the APTAMD refined model shows better agreement with the NMR conformation of the central segment than the initial RNAComposer structure.     

<img src="./Images/aptamd_label.png" width="450" height="225" style="display: block; margin: 0 auto">
</div>


### Other script tools in APTAMD

The APTAMD suite includes scripts for MM parameterization,  system preparation and edition, entropy calculations, Autodock calculations on multiple MD frames of receptors and ligands, etc. These and other scripts will be uploaded to GitHub (and documented) in the near future. 


```python

```
