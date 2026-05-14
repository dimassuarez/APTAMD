
Maximum numer of atoms is now 25000

Files bresna.F90 ,  nmode.F90 , eprnt.F90 and  forces.F90 are modified for belly NMODE.
Look for comment lines headed by  ! Di+

Big integers (-i8) are needed to treat large systems.

When comparing NMODE and SANDER outputs with BELLY-fixed waters, 
BOND, ANGLE, DIHED, NB14 and  ELEC14 terms are identical. 
There is a difference in the NB and ELEC terms beacuse NMODE Belly
implementation does not delete non-bonding interactions among 
the BELLY waters. Nonetheless, everything is OK as long as 
we minimize with SANDER TCNG and then do the Hessian calc with NMODE.  

RMS grad value is very sensitive to the XYZ coordinates. Thus, SANDER optimized
geometry must be saved in NETCDF format. Then this NETCDF file is transformed
into an ascii file (XYZ coords have thus a lot of digits). The modified  NMODE
can read this XYZ file (ntx=2)

Other changes have been made to allow a variable 1-4 scaling (needed for GLYCAM).
See rdparm.f, alloc.F90, Ephi.F  and forces.F90 


