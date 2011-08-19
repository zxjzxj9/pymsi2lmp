#!/usr/bin/python
import re
from numpy import abs

""" 
 frc.py - this module provides a method to update an frc file with new
   parameters, where the new parameters are concatenated in a list.
"""

# Compass header lines matched up with missing types from converter.py.
headers = """
#quartic_bond           b
#quartic_angle          a
#bond-bond              bb
#bond-bond_1_3          bb13
#bond-angle             ba
#torsion_3              tor
#end_bond-torsion_3     ebt
#middle_bond-torsion_3  mbt
#angle-torsion_3        at
#wilson_out_of_plane    oop
#angle-angle            aa
#angle-angle-torsion_1  aat
#nonbond(9-6)           vdw
"""
## Puts headers in a dictionary that maps to their shortened type.
headers = dict([x.split() for x in headers.strip().split('\n')])

## Should equilibrium torsion angles be allowed to be nonzero.
nonzero_torsion_angles = False
nonzero_oop_angles = False


def update(oldfrc, newfrc, unknown, p):
    original = open(oldfrc, 'r')
    modified = open(newfrc, 'w')
    
    # Regex that will match only a data line.
    # Looks for a line with "#.#  #"
    data_re = re.compile('^\s*\d\.\d\s*\d.+\n$')
    ct = 0
    while True:
        # Read until we encounter the end of the file.
        line = original.readline()
        if line == '': break
        modified.write(line)
        cols = line.split()
        if line[0] == '#' and cols[0] in headers:
            # Wait until we have hit the first actual data point.
            while True:
                line = original.readline()
                if not data_re.match(line):
                    modified.write(line)
                    continue
                # First data point found, now write new terms.
                term = headers[cols[0]]
                if term in unknown: 
                    for atom_group in unknown[term]:
                        ct = write_atom_group(modified, term, atom_group, p, ct)
                modified.write(line)
                break

# Writes a formatted parameter line to the frc file.
def write_atom_group(fid, term, atoms, p, ct):
    if term == 'a' and abs(p[ct])<10.0 : p[ct] = 120.0
    if term == 'b' and p[ct]<0.1: p[ct] = 1.5
    natoms = len(atoms)
    line = ' 1.0   0   ' + natoms*('  %-4s') + ' '
    line = line%tuple(atoms)
    n =  term_parameter_count(term, atoms)
    if term == 'tor' and not nonzero_torsion_angles:        
        line += n*'   %9.4f    0.0'       
    elif term == 'oop' and not nonzero_oop_angles:
        line += n*'    %10.4f    0.0'       
    else:
        line += n*'   % 9.4f'
       
    line = line%tuple(p[ct:ct+n])
    ct += n
    fid.write(line + '\n')    
    print '%4s'%term, line
    return ct

# Returns parameters count for a term, some require atoms due to symmetry.
def term_parameter_count(term, atoms=[]):    
    sizes = {'a':4,'bb':1,'aat':1,'tor':3,
             'b':4,'aa':1,'mbt':3,'oop':1, 'bb13':1}
    if nonzero_torsion_angles: sizes['tor'] = 6
    if nonzero_oop_angles:     sizes['oop'] = 2
    if term in sizes: return sizes[term]
    # Bond angle can be unsymmetric if i!=k for ijk.
    elif term == 'ba':
        if atoms[0] == atoms[2]: return 1
        return 2
    # Torsions have double constants if not mirror symmetric. 
    elif term == 'ebt':
        if atoms[0] != atoms[3] or atoms[1] != atoms[2]: return 6
        return 3
    # Angle torsion has double constants if not mirror symmetric
    elif term == 'at':
        if atoms[0] != atoms[3] or atoms[1] != atoms[2]: return 6
        return 3
    else: print 'Unknown parameter', term, atoms
