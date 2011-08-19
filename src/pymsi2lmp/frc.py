#!/usr/bin/env python
"""
 frc.py - this module provides methods to read a COMPASS frc file and provide
          parameters for an interaction for a specific set of atom types.
          Interactions are alphabetically sorted by atom type in such that
          their symmetries are preserved, e.g. dihedral DCBA becomes ABCD.
          Lists of atom types are stored by ':' delimited strings so that they
          are hashable and can be stored in a dictionary for fast searching.          
"""

# Abbreviations for each of the types of interactions in the COMPASS potential.
# Each interaction is of class: NONBOND, BOND, ANGLE, TORSION, or OOP.
compass_key = { 'b':    ['#quartic_bond',          'BOND'],
                'a':    ['#quartic_angle',         'ANGLE'],
                'bb':   ['#bond-bond',             'ANGLE'],
                'bb13': ['#bond-bond_1_3',         'TORSION'],
                'ba':   ['#bond-angle',            'ANGLE'],
                'tor':  ['#torsion_3',             'TORSION'],
                'ebt':  ['#end_bond-torsion_3',    'TORSION'],
                'mbt':  ['#middle_bond-torsion_3', 'TORSION'],
                'at':   ['#angle-torsion_3',       'TORSION'],
                'oop':  ['#wilson_out_of_plane',   'OOP'],
                'aa':   ['#angle-angle',           'OOP'],
                'aat':  ['#angle-angle-torsion_1', 'TORSION'],
                'vdw':  ['#nonbond(9-6)',          'NONBOND']}

# Returns a file handle at the beginning of a section in an frc file.
def open_at(path, section):
    try:
        fid = open(path, 'r')  
        for line in fid: 
            if line.startswith(section): return fid    
    except:
        print 'Error: frc file,', path, 'cannot be opened.'
        import sys
        sys.exit()
    
    print 'Error: frc file does not have', section
    import sys
    sys.exit()
                
# Returns the sorting function used by the interaction.
def sort_function(interaction):
    style = compass_key[interaction][1]
    # Special case, AA sorts only forward/back.
    if interaction =='aa': return sort_aa
    if   style=='BOND' or style=='ANGLE': return sort_bond_or_angle
    elif style=='OOP':                    return sort_oop
    elif style=='TORSION':                return sort_torsion
    else:                                 return lambda x:x    
    
# Bonds and angles are reversible AB = BA, ABC = CBA.
def sort_bond_or_angle(types, indices=None):
    if types[0] > types[-1]:
        types = types[::-1]
        if indices: indices = indices[::-1]
    if indices: return types, indices
    return types

# Dihedrals are symmetric ABCD = DCBA.
def sort_torsion(types, indices=None):
    # If first and last are the same, then check middle two atoms.
    if types[0]>types[3] or (types[0]==types[3] and types[1]>types[2]):
        types = types[::-1]
        if indices: indices = indices[::-1]

    if indices: return types,indices
    return types

# Impropers are symmetric about B in ABCD.
def sort_oop(types, indices=None):
    def swap_types_and_indices(i, j, types, indices):
        if types[i] < types[j]: 
            types[i],types[j] = types[j],types[i]
            if indices:
                indices[i],indices[j] = indices[j],indices[i]     
        return types, indices

    types,indices = swap_types_and_indices(0, 2, types, indices)
    types,indices = swap_types_and_indices(0, 3, types, indices)
    types,indices = swap_types_and_indices(2, 3, types, indices)
    if indices: return types,indices
    return types

# Angle-angle can swap first and last.
def sort_aa(types, indices=None):
    if types[0] > types[3]:
        types[0],types[3] = types[3],types[0]
        if indices:
            indices[0],indices[3] = indices[3],indices[0]
    if indices: return types, indices
    return types
    
def next_section(s): return s.startswith('#') and s.endswith('compass\n')
def skip(s): return s=='\n' or s[0]=='!' or s[0]=='@' or s[0]=='>'

# Reads an frc file and supplies the parameters needed for a lammps input.
class Frc:
    def __init__(self, path):
        self.path  = path
        self.types = self.read_masses()
        self.read_equiv()
                
        self.coeff = {}        
        for i in compass_key:
            self.coeff[i] = self.readparam(i)
  
    # Finds forcefield coefficients from the interaction type and 
    # forcefield types of a set of atoms.
    def get_param(self, fftypes, interaction):                        
        ffstr = ':'.join(fftypes)
        coeff = self.coeff[interaction]                                
        if ffstr in coeff:
            return coeff[ffstr]
        
        sortf = sort_function(interaction)
        # Try again with equivalent types.  NOTE: tries to replace all atoms.
        if   compass_key[interaction][1]=='NONBOND': style = 1           
        elif compass_key[interaction][1]=='BOND':    style = 2
        elif compass_key[interaction][1]=='ANGLE':   style = 3
        elif compass_key[interaction][1]=='TORSION': style = 4
        elif compass_key[interaction][1]=='OOP':     style = 5
        
        
        eqffstr = ':'.join(sortf([self.equiv[a][style-1] for a in fftypes]))

        # Optional: if you want to see which equivalences are being used.
        # This can end up producing a lot of output.
        #print 'Replacing', ffstr, 'with', eqffstr, 'for', interaction        
        if eqffstr in coeff: return coeff[eqffstr]
        
        if compass_key[interaction][1]=='TORSION':
            # Ok - now we are desperate - try wildcards.
            ffstr = ':'.join(sortf(fftypes[:-1]+['*']))        
            if ffstr in coeff: return coeff[ffstr]
            ffstr = ':'.join(sortf(['*']+fftypes[1::]))
            if ffstr in coeff: return coeff[ffstr]
            ffstr = ':'.join(sortf(['*']+fftypes[1:-1]+['*']))
            if ffstr in coeff: return coeff[ffstr]
        return None

    # Reads the equivalences table.
    def read_equiv(self):
        fid = open_at(self.path, '#equivalence')
        self.equiv = {}
        for line in fid:
            if skip(line):         continue
            if next_section(line): break
            x = line.split()
            self.equiv[x[2]] = x[3::]    
        
    # Reads the atomic masses field.
    def read_masses(self):
        fid = open_at(self.path, '#atom_types')
        types = {}        
        for line in fid:
            if skip(line):         continue
            if next_section(line): break
            # Remove comment from line and split into columns.
            x = line[0:39].split() 
            if len(x) != 5:        continue
            # types[fftype] = element, mass.
            types[x[2]] = [x[4], float(x[3])] 
        return types
            
    # Reads parameters from each interaction type.
    def readparam(self, interaction):
        table = {}        
        fid   = open_at(self.path, compass_key[interaction][0])
        sortf = sort_function(interaction)
        for line in fid:
            if skip(line):         continue
            if next_section(line): break
            x = line.split()
            # Determines how many atoms types are listed in the line.
            # All fftypes start with a letter or * for wildcard.
            # Coefficients will not (as they should start with a number or -).
            # One could also have a lookup table by interaction type.
            for i in range(1,5):
                if not (x[2+i][0].isalpha() or x[2+i][0]=='*'):
                    natoms = i
                    break
            fftypes = x[2:2+natoms]
            param   = [float(y) for y in x[2+natoms::]]
            s_fftypes = sortf(fftypes[:])
            
            # Some interactions coefficients are not symmetric.
            # Here we make sure they are not messed up with the sorting.
            # If bond angle was reversed, then reverse parameters.
            if interaction == 'ba'  and fftypes[0] != s_fftypes[0]:
                param = param[::-1]                                                        
            # If end-bond-torsion or angle-torsion was reversed, 
            # swap left and right parameters.
            if any(interaction==x for x in ['ebt','at']) and fftypes[0]!=s_fftypes[0]:
                param[0:3],param[3:6] = param[3:6],param[0:3]                                                
            table[':'.join(s_fftypes)] = param
        return table
