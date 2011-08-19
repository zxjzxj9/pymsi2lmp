#!/usr/bin/env python
"""
 pymsi2lmp.py - main function for converting models from Materials Studio
                into a LAMMPS readable format.
    AUTHOR
        Jay Oswald: j-oswald@asu.edu
                
    SYNOPSIS
        pymsi2lmp.py [-i INPUTPATH] [-frc COMPASSPATH]
        
    DESCRIPTION
        Converts INPUTPATH.mdf and INPUTPATH.car to INPUTPATH.lammps using the
        COMPASS forcefield parameters specified by COMPASSPATH.
        
        -i INPUTFILE
            Basename of files exported by Materials Studio.  If left blank, the
            first *.mdf found by glob will be used.
            
        -frc COMPASS PATH
            Specifies the location of the frc parameter file.  If compass.frc is
            located in the working path, then it will be used as a default.    
"""
import glob
import insight
import lammps_writer
import frc

def msi2lmp(rootname, frcpath):
    system  = insight.get_system(rootname) 
    compass = frc.Frc(frcpath)
    missing = lammps_writer.write_data(system, compass)
    return missing

def main(args):
    # Default frc path.
    frcpath = 'compass.frc'
    # Sets a default rootname.
    mdf = glob.glob('*.mdf')
    if len(mdf) > 0: rootname = mdf[0][:-4]

    # Input arguments
    for i in range(1, len(args)):
        if args[i] == '-i':     rootname = args[i+1]
        elif args[i] == '-frc': frcpath  = args[i+1]
        else: continue

    missing = msi2lmp(rootname, frcpath)
    for m in missing: 
        term = frc.compass_key[m[0]][0][1::]
        print 'Unable to find', term, 'data for', ' '.join(m[1])
    
    
if __name__ == '__main__':
    import sys
    main(sys.argv)
