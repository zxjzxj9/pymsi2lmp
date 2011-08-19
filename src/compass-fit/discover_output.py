#!/usr/bin/env python
"""
    Reads output files created by MS discover.
"""

from numpy import array

# Reads atom coords and energy in the materials studio trajectory file.
# Each timestep trajectory is packed into a 1D array
# with the dimension varying fastest, eg [x1 y1 z1 ... xn yn zn].
def read_arc_file(base):
    arc = open(base+ '.arc', 'r')
    energy     = []
    trajectory = []
    while True:
        line = arc.readline().strip()
        if not line: break
        if line.startswith('Materials Studio Generated CAR File'):
            energy.append(float(line.split()[-1]))
            trajectory.append([])
            # Ignore !DATE line
            arc.readline()
            while True:
                line = arc.readline().strip()
                if not line or line=='end': break
                for s in line.split()[1:4]:
                    trajectory[-1].append(float(s))
    return trajectory, array(energy)

# Reads the MS disco energy file (useful if you want to try to compare term-by-term.
def read_disco_energy_file(base):    
    e = {'ba':0.0,'bb':0.0,'mbt':0.0,'ebt':0.0,'aat':0.0,'aa':0.0,'oop':0.0,
         'bb13':0.0,'at':0.0}
    for line in open(base + '.out', 'r'):
        cols = line.split()
        if len(cols)!=2: pass
        elif cols[0] == 'Total:':               e['pe']   = float(cols[1])
        elif cols[0] == 'Bond:':                e['b']    = float(cols[1])
        elif cols[0] == 'Angle:':               e['a']    = float(cols[1])
        elif cols[0] == 'Torsion:':             e['tor']  = float(cols[1])
        elif cols[0] == 'OutOfPlane:':          e['oop']  = float(cols[1])
        elif cols[0] == 'BondBond:':            e['bb']   = float(cols[1])
        elif cols[0] == 'BondAngle:':           e['ba']   = float(cols[1])
        elif cols[0] == 'EndBondTorsion:':      e['ebt']  = float(cols[1])
        elif cols[0] == 'MiddleBondTorsion:':   e['mbt']  = float(cols[1])
        elif cols[0] == 'AngleTorsion:':        e['at']   = float(cols[1])
        elif cols[0] == 'AngleAngleTorsion:':   e['aat']  = float(cols[1])
        elif cols[0] == 'BondBond-1-3:':        e['bb13'] = float(cols[1])
        elif cols[0] == 'AngleAngle:':          e['aa']   = float(cols[1])
        elif cols[0] == 'Vdw:':                 e['vdw']  = float(cols[1])
        elif cols[0] == 'Electrostatic:':       e['q']    = float(cols[1])
    return e        