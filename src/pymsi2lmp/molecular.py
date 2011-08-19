#!/usr/bin/env python
"""
 molecular.py: this module provides the basic containers needed 
               to build lammps data files.
"""

from frc import sort_bond_or_angle, sort_oop, sort_torsion

# 
class AtomSet:
    def __init__(self, atoms, type_index):
        self.atoms      = atoms
        self.type_index = type_index        

# Atomic data class.
class Atom:
    def __init__(self):
        self.x    = []
        self.conn = []
        self.seq  = 0
        self.ff   = 0
        self.sym  = 0        
    # Writes to a string.
    def __str__(self):
        s = '%5s: % 9.6f % 9.6f % 9.6f' %tuple([self.ff] + self.x)
        s += '  ' + str(self.conn)
        return s

# System of atoms and bonds.
class System:
    def __init__(self, atoms, title, pbc, bounds):
        self.atoms  = atoms     # List of atoms in the system.
        self.title  = title     # System title.        
        self.pbc    = pbc       # Periodic boundary conditions.
        if bounds==None: self.bounds = self.compute_bounds()
        else:            self.bounds = bounds

    # Remaps atom coordinates so that they fit in the simulation box.
    def remap_to_box(self):
        if not self.pbc: return
        dx = [self.bounds[1]-self.bounds[0], 
              self.bounds[3]-self.bounds[2], 
              self.bounds[5]-self.bounds[4]]
        for a in self.atoms:
            for i in range(3):
                while a.x[i] < self.bounds[i*2]:
                    a.x[i] += dx[i]
                while a.x[i] > self.bounds[i*2+1]:
                    a.x[i] -= dx[i]
                            
    # Computes the bounds of the atoms in the domain. 
    def compute_bounds(self):
        r = []
        for i in range(3):
            x  = [a.x[i] for a in self.atoms]
            r += [min(x), max(x)]
        return tuple(r)

    # Returns an array of atomic types.
    def atom_types(self):
        types = []
        for a in self.atoms:
            if not a.ff in types:
                types.append(a.ff)
        return self.atoms,types

    # Builds the bond table up.
    def bonds(self):
        types,bonds = [],[]
        for i, a in enumerate(self.atoms):
            for j in a.conn:
                if i > j: continue
                ij = [i,j]
                t, ij = sort_bond_or_angle([self.atoms[p].ff for p in ij], ij)
                if not t in types: types.append(t)
                bonds.append(AtomSet(ij, types.index(t)))
        return bonds, types
    
    # Builds up the angles from each atom.
    def angles(self):
        types,angles = [],[]
        # Center is atom j.
        for j, a in enumerate(self.atoms):
            for i in a.conn:
                for k in a.conn:
                    if k <= i: continue
                    ijk = [i,j,k]
                    t,ijk = sort_bond_or_angle([self.atoms[p].ff for p in ijk], ijk)
                    if not t in types: 
                        types.append(t)
                    angles.append(AtomSet(ijk, types.index(t)))
        return angles, types

    # Builds up improper groups.
    def impropers(self):
        types,oop = [],[]
        # Center is atom j
        for j, a in enumerate(self.atoms):
            if len(a.conn)==3:
                ijkl = [a.conn[0], j, a.conn[1], a.conn[2]]
                t,ijkl = sort_oop([self.atoms[p].ff for p in ijkl], ijkl)
                if not t in types: 
                    types.append(t)
                oop.append(AtomSet(ijkl, types.index(t)))
        # Now loop over sets of 4 atoms.
        for j, a in enumerate(self.atoms):        
            if len(a.conn)==4:                
                all_oop = [[a.conn[1], j, a.conn[2], a.conn[3]],
                           [a.conn[0], j, a.conn[2], a.conn[3]],
                           [a.conn[0], j, a.conn[1], a.conn[3]],
                           [a.conn[0], j, a.conn[1], a.conn[2]]]                
                all_oop = [sort_oop([self.atoms[p].ff for p in q], q) for q in all_oop]
                for t,ijkl in all_oop:
                    if not t in types:
                        types.append(t)
                    oop.append(AtomSet(ijkl, types.index(t)))
        return oop, types

    # Builds up dihedral groups.
    def dihedrals(self):
        types,dihed = [],[]
        for j, a in enumerate(self.atoms):
            if len(a.conn)<2: continue
            for k in a.conn:
                if not j<k or len(self.atoms[k].conn)<2: continue
                for i in self.atoms[j].conn:
                    if i==k: continue
                    for l in self.atoms[k].conn:
                        if l==j: continue
                        ijkl = [i,j,k,l]
                        t,ijkl = sort_torsion([self.atoms[p].ff for p in ijkl], ijkl)
                        if not t in types:
                            types.append(t)
                        dihed.append(AtomSet(ijkl, types.index(t)))
        return dihed, types
        