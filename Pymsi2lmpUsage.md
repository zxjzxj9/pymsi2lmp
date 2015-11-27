# pymsi2lmp instructions #
This is a drop in replacement for the c version of msi2lmp provided with LAMMPS.

# Usage #
**pymsi2lmp.py `[`-i INPUTPATH`]` `[`-frc COMPASSPATH`]`**

## DESCRIPTION ##
> Converts INPUTPATH.mdf and INPUTPATH.car to INPUTPATH.lammps using the COMPASS forcefield parameters specified by COMPASSPATH.

> -i INPUTFILE
> > Basename of files exported by Materials Studio.  If left blank, the first `*`.mdf found by glob will be used.


> -frc COMPASS PATH
> > Specifies the location of the frc parameter file.  If compass.frc is
> > located in the working path, then it will be used as a default.