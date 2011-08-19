#!/usr/bin/env python
"""
 lammps_writer.py - These routines export the system and a set of forcefield 
                    parameters from a frc file into an input file for LAMMPS
                    with a CLASS2 potential.  Potentially, different types of 
                    frc files (e.g. pcff) could be used.
"""

import frc2lmp

# Template string for the system bounds specified to LAMMPS.
lmp_bounds = """
 %15.9f %15.9f xlo xhi
 %15.9f %15.9f ylo yhi
 %15.9f %15.9f zlo zhi
"""

# Writes the LAMMPS output file. 
def write_data(system, frc):
    fid = open(system.title+'.lammps','w')
    atoms,  types  = system.atom_types()
    bonds,  btypes = system.bonds()
    angles, atypes = system.angles()
    dihed,  dtypes = system.dihedrals()
    oop,    otypes = system.impropers()       
        
    fid.write('LAMMPS 2005 data file for ' + system.title + '\n\n')
    fid.write(' %6d atoms\n'       %len(atoms))
    fid.write(' %6d bonds\n'       %len(bonds))
    fid.write(' %6d angles\n'      %len(angles))
    fid.write(' %6d dihedrals\n'   %len(dihed))
    fid.write(' %6d impropers\n\n' %len(oop))         
    fid.write(' %3d atom types\n'  %len(types))
    fid.write(' %3d bond types\n'  %len(btypes))
    fid.write(' %3d angle types\n' %len(atypes))          
    if len(dtypes): fid.write(' %3d dihedral types\n' %len(dtypes))
    if len(otypes): fid.write(' %3d improper types\n' %len(otypes))
    fid.write(lmp_bounds %system.compute_bounds())
    
    fid.write('\nMasses\n\n')
    for i,t in enumerate(types):
        fid.write('%4d %10.6f\n' %(i+1, frc.types[t][1]))

    fid.write('\nPair Coeffs\n\n')       
    vdwcoeff,vdwmissing = frc2lmp.pair(types, frc)
    for i,c in enumerate(vdwcoeff):
        fid.write(' %3d'%(i+1) + len(c)*' %14.10f'%c + '\n')
    
    bondcoeff,bondmissing = frc2lmp.bond(btypes, frc)
    if len(btypes): fid.write('\nBond Coeffs\n\n')
    for i,c in enumerate(bondcoeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')
       
    anglecoeff,anglemissing = frc2lmp.angle(atypes, frc)
    if len(atypes): fid.write('\nAngle Coeffs\n\n')
    for i,c in enumerate(anglecoeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')        
   
    torsioncoeff,torsionmissing = frc2lmp.torsion(dtypes, frc)
    if len(dtypes): fid.write('\nDihedral Coeffs\n\n')
    for i,c in enumerate(torsioncoeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')
    
    oopcoeff,oopmissing = frc2lmp.oop(otypes, frc)
    if len(otypes): fid.write('\nImproper Coeffs\n\n')
    for i,c in enumerate(oopcoeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')

    bbcoeff,bbmissing = frc2lmp.bondbond(atypes,frc)
    if len(atypes): fid.write('\nBondBond Coeffs\n\n')
    for i,c in enumerate(bbcoeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')
        
    bacoeff,bamissing = frc2lmp.bondangle(atypes,frc)
    if len(atypes): fid.write('\nBondAngle Coeffs\n\n')
    for i,c in enumerate(bacoeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')
        
    aacoeff,aamissing = frc2lmp.angleangle(otypes,frc)
    if len(otypes): fid.write('\nAngleAngle Coeffs\n\n')
    for i,c in enumerate(aacoeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')
    
    aatcoeff,aatmissing = frc2lmp.angleangletorsion(dtypes,frc)
    if len(dtypes): fid.write('\nAngleAngleTorsion Coeffs\n\n')
    for i,c in enumerate(aatcoeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')

    ebtcoeff,ebtmissing = frc2lmp.endbondtorsion(dtypes,frc)
    if len(dtypes): fid.write('\nEndBondTorsion Coeffs\n\n')
    for i,c in enumerate(ebtcoeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')
    
    mbtcoeff,mbtmissing = frc2lmp.midbondtorsion(dtypes,frc)
    if len(dtypes): fid.write('\nMiddleBondTorsion Coeffs\n\n')
    for i,c in enumerate(mbtcoeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')
        
    bb13coeff,bb13missing = frc2lmp.bondbond13(dtypes,frc)
    if len(dtypes): fid.write('\nBondBond13 Coeffs\n\n')
    for i,c in enumerate(bb13coeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')
        
    atcoeff,atmissing = frc2lmp.angletorsion(dtypes,frc)
    if len(dtypes): fid.write('\nAngleTorsion Coeffs\n\n')
    for i,c in enumerate(atcoeff):
        fid.write('%3d'%(i+1) + len(c)*' %10.4f'%c + '\n')
        
    missing  = vdwmissing + bondmissing + anglemissing + torsionmissing
    missing += oopmissing + bbmissing   + bamissing    + aamissing + aatmissing
    missing += ebtmissing + mbtmissing  + bb13missing  + atmissing       
    
    #  Now list atoms, bonds, angles, etc.
    fid.write('\nAtoms\n\n')
    for i,x in enumerate(atoms):
        fid.write(' %6d %6d %3d %9.6f' %(i+1, x.seq, types.index(x.ff)+1, x.q))
        fid.write('%15.9f %15.9f %15.9f\n' %tuple(x.x))

    # Writes a string of atom indices to a file.
    def write_atom_indices(fid, atoms):
        fid.write(('%6d '*len(atoms)).strip() %tuple([a+1 for a in atoms]) + '\n')
        
    fid.write('\nBonds\n\n')
    for i,x in enumerate(bonds):
        fid.write('%6d %3d' %(i+1, x.type_index+1))
        write_atom_indices(fid, x.atoms)

    fid.write('\nAngles\n\n')
    for i,x in enumerate(angles):
        fid.write('%6d %3d ' %(i+1, x.type_index+1))
        write_atom_indices(fid, x.atoms)

    fid.write('\nDihedrals\n\n')
    for i,x in enumerate(dihed):
        fid.write('%6d %3d ' %(i+1, x.type_index+1)) 
        write_atom_indices(fid, x.atoms)
        
    fid.write('\nImpropers\n\n')
    for i,x in enumerate(oop):
        fid.write('%6d %3d ' %(i+1, x.type_index+1)) 
        write_atom_indices(fid, x.atoms)
    
    return missing