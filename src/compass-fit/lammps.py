#!/usr/bin/env python
"""
 This script makes a LAMMPS input file to fit a compass potential.
"""

import subprocess
from subprocess import PIPE
default_lmp_input = """echo none
units               real
boundary            s s s 
atom_style          full
pair_style          lj/class2/coul/cut 9.5 9.5
%s
pair_modify         shift no 
special_bonds       lj/coul 0.0 0.0 1.0 angle yes dihedral yes
read_data           %s
#neighbor            0.3 bin
neighbor            0.3 nsq
thermo_style        custom pe evdwl ecoul ebond eangle edihed eimp
fix                 1 all nvt temp 300 300 1
timestep            0
run                 1
"""

# Scans the input file to see if pair/bond/angle/dihedral/improper are there.
def input_file_terms(inputfile):
    types = '' 
    for line in open(inputfile, 'r'):
        line = line.strip()
        if   line.endswith('bond types'):     types += 'bond_style  class2\n'
        elif line.endswith('angle types') and int(line[0]):  
            types += 'angle_style class2\n'
        elif line.endswith('dihedral types'): types += 'dihedral_style class2\n'
        elif line.endswith('improper types'): types += 'improper_style class2\n'
        elif line.endswith('xlo xhi'): return types
    return types

# Writes the lammps input file.
def run(base, step, lmp):
    inputfile = base + '.lammps.%03d'%step
    p = subprocess.Popen([lmp], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    types = input_file_terms(inputfile)
    return p.communicate(default_lmp_input %(types, inputfile))[0]

# Extracts the potential energy computed by LAMMPS.
def extract_energy(out):
    beg = out.find('PotEng')    
    end = out.find('Loop time', beg)-1
    if end < 0: 
        import sys
        print out
        sys.exit(1)
    out = out[beg:end]
    beg = out.rfind('\n')+1
    energy_line = out[beg:end] 
    
    return float(energy_line.split()[0])

# Reads the number of atoms from a lammps data file.
def get_num_atoms(base):
    tag = ' atoms'
    for line in open(base + '.lammps05'):
        line = line.strip()
        if line.endswith(tag):
            return int(line[0:-len(tag)])

# Reads the data file and modifies the atomic positions.
# Each new file is saved as base.lammps.05.step
def modify_data_file(base, x, step):
    old = open(base + '.lammps', 'r')
    new = open(base + '.lammps.%03d'%step, 'w')

    def bounds(d, lo, hi):
        return ' %15.9f %15.9f %slo %shi\n' %(lo,hi,d,d)

    # Copy line for line until encounter the Atoms line.
    # Replace the bounds with the updated bounds.
    for line in old:
        if line.strip().endswith('xlo xhi'):
            new.write(bounds('x', min(x[0::3]), max(x[0::3])))
        elif line.strip().endswith('ylo yhi'):
            new.write(bounds('y', min(x[1::3]), max(x[1::3])))
        elif line.strip().endswith('zlo zhi'):
            new.write(bounds('z', min(x[2::3]), max(x[2::3]))) 
        else: new.write(line)
        if line.startswith('Atoms'): break 

    # Write new atom lines (or add blank lines).
    i = 0
    for line in old:
        if line.strip() == '': 
            new.write(line)
            continue
        elif line.startswith('Bonds'):
            new.write(line)
            break
        c = line.split()        
        c = (int(c[0]),int(c[1]),int(c[2]),float(c[3]),x[i],x[i+1],x[i+2])
        new.write(' %6d %6d %3d %9.6f %15.9f %15.9f %15.9f\n' %c)
        i += 3
    # Copy remainder of file.
    for line in old: new.write(line)

# Reads the potential energies reported in the log file.
def get_pe_from_log_file(log=''):
    if log == '': log = 'log.lammps'
    fid = open(log, 'r')
    while True:
        line = fid.readline()
        if line == '': break
        if line.startswith('PotEng'): 
            variables = line.split()
            energies  = [float(x) for x in fid.readline().split()]              
    return dict(zip(variables, energies))
    
    
