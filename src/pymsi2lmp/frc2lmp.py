#!/usr/bin/env python
"""
 frc2lmp.py - These routines are for converting frc parameter data into a 
              format that can easily be written to a LAMMPS input file.
"""
from frc import sort_bond_or_angle as sort
from frc import sort_aa as sortaa
# Returns the pairwise coefficients in an list for each type.
def pair(types, frc):
    coeffs,missing = [],[]
    for t in types:
        if t in frc.coeff['vdw']:
            c = frc.coeff['vdw'][t][::-1]
        else:
            missing += [['vdw', t]]
            c = 2*[0.0]
        coeffs.append(tuple(c))
    return coeffs,missing
        
# Returns the bond coefficients in an list for each type.       
def bond(types, frc):
    coeffs,missing = [],[]
    for t in types:
        c = frc.get_param(t, 'b') 
        if c == None:
            missing += [['b', t]]
            c = 4*[0.0]        
        coeffs.append(tuple(c))
    return coeffs,missing

# Returns the angle coefficients in an list for each type.
def angle(types, frc):
    coeffs,missing = [],[]
    for t in types:        
        c = frc.get_param(t, 'a')        
        if c == None:
            missing += [['a', t]]
            c = 4*[0.0]
        coeffs.append(tuple(c))
    return coeffs,missing
        
# Returns the torsion coefficients in an list for each type.
def torsion(types, frc):
    coeffs,missing = [],[]
    for t in types:
        c = frc.get_param(t, 'tor')
        if c == None:
            missing += [['tor', t]]
            c = 6*[0.0]
        coeffs.append(tuple(c))
    return coeffs,missing

# Returns the  coefficients in an list for each type.
def oop(types, frc):
    coeffs,missing = [],[]
    for t in types:
        c = frc.get_param(t, 'oop')
        if c == None:
            missing += [['oop', t]]
            c = 2*[0.0]
        coeffs.append(tuple(c))
    return coeffs,missing

# Returns the bond-bond coefficients in an list for each type.
def bondbond(types, frc):
    coeffs,missing = [],[]
    for t in types:        
        b1,b2 = frc.get_param(sort(t[0:2]),'b'),frc.get_param(sort(t[1:3]),'b')        
        bb = frc.get_param(t, 'bb')
        try:
            c = bb + [b1[0], b2[0]]
        except:
            missing += [['bb',t]]
            c = 3*[0.0]        
        coeffs.append(tuple(c))
    return coeffs,missing
        
# Returns the bond-angle coefficients in an list for each type.
def bondangle(types, frc):
    coeffs,missing = [],[]
    for t in types:        
        b1,b2 = frc.get_param(sort(t[0:2]),'b'),frc.get_param(sort(t[1:3]),'b')
        ba = frc.get_param(t, 'ba')
        try:
            if len(ba) == 1: ba += ba
            c = ba + [b1[0], b2[0]]
        except:                            
            c = 4*[0.0]        
            missing += [['ba',  t]]
        coeffs.append(tuple(c))
    return coeffs,missing
        
# Returns the angle-angle coefficients in an list for each type.
def angleangle(types, frc):
    coeffs,missing = [],[]        
    for t in types:        
        #   jj'_______k'
        #    /\
        #   /  \        E = K*(theta-theta0)*(theta'-theta0')
        # i/    \k,i'   
        #                ijkl                  kjil                        ijlk
        aa1,aa2,aa3 = sortaa(t[:]),sortaa([t[3],t[1],t[0],t[2]]),sortaa([t[0],t[1],t[3],t[2]])                
        m1,m2,m3 = frc.get_param(aa1,'aa'),frc.get_param(aa2,'aa'),frc.get_param(aa3,'aa')       
        
        # Gives the three types of angles possible in a torsion term.
        #            a1=ijk                a2=ijl                  a3=kjl        
        a1,a2,a3 = sort(t[0:3]), sort([t[0],t[1],t[3]]), sort([t[2],t[1],t[3]])
        t1,t2,t3 = frc.get_param(a1,'a'),frc.get_param(a2,'a'),frc.get_param(a3,'a')        
        
        if any([p==None for p in [t1,t2,t3,m1,m2,m3]]):                    
            if m1==None and not ['aa',aa1] in missing: missing += [['aa',aa1]]
            if m2==None and not ['aa',aa2] in missing: missing += [['aa',aa2]]
            if m3==None and not ['aa',aa3] in missing: missing += [['aa',aa3]]
            
            if t1 == None: missing += [[ 'a', a1]]
            if t2 == None: missing += [[ 'a', a2]]
            if t3 == None: missing += [[ 'a', a3]]
            c = 6*[0.0]
        else:
            c = m1 + m2 + m3 + [t1[0],t2[0],t3[0]]
        coeffs.append(tuple(c))

    return coeffs,missing
    
# Returns the angle-angle-torsion coefficients in an list for each type.
def angleangletorsion(types, frc):
    coeffs,missing = [],[]
    for t in types:
        a1,a2=sort(t[0:3]),sort(t[1:4])
        t1,t2 = frc.get_param(a1,'a'),frc.get_param(a2,'a')
        aat   = frc.get_param(t, 'aat')                
        if t1==None or t2==None or aat==None:
            missing += [['aat', t]]
            c = 3*[0.0]
        else:
            c = aat + [t1[0], t2[0]]        
        coeffs.append(tuple(c))
    return coeffs,missing

# Returns the end-bond-torsion coefficients in an list for each type.
def endbondtorsion(types, frc):
    coeffs,missing = [],[]
    for t in types:
        b1,b3 = frc.get_param(sort(t[0:2]),'b'), frc.get_param(sort(t[2:4]),'b')
        ebt   = frc.get_param(t, 'ebt')
        if b1==None or b3==None or ebt==None:
            missing += [['ebt', t]]
            c = 8*[0.0]
        else:
            if len(ebt) == 3: ebt += ebt
            c = ebt + [b1[0], b3[0]]
        coeffs.append(tuple(c))
    return coeffs,missing
    
# Returns the mid-bond-torsion coefficients in an list for each type.
def midbondtorsion(types, frc):
    coeffs,missing = [],[]
    for t in types:
        mb  = frc.get_param(sort(t[1:3]),'b')
        mbt = frc.get_param(t, 'mbt')
        if mb==None or mbt==None:
            missing += [['mbt', t]]
            c = 4*[0.0]
        else:            
            c = mbt + [mb[0]]
        coeffs.append(tuple(c))
    return coeffs,missing

# Returns the bond-bond13 coefficients in an list for each type.    
def bondbond13(types, frc):
    coeffs,missing = [],[]
    for t in types:
        b1,b3 = frc.get_param(sort(t[0:2]),'b'), frc.get_param(sort(t[2:4]),'b')        
        bb13  = frc.get_param(t, 'bb13')
        if b1==None or b3==None or bb13==None:
            missing += [['bb13', t]]
            c = 3*[0.0]
        else:            
            c = bb13 + [b1[0], b3[0]]
        coeffs.append(tuple(c))
    return coeffs,missing
    
# Returns the mid-bond-torsion coefficients in an list for each type.
def angletorsion(types, frc):
    coeffs,missing = [],[]
    for t in types:
        a1,a2 = frc.get_param(sort(t[0:3]),'a'),frc.get_param(sort(t[1:4]),'a')        
        at    = frc.get_param(t, 'at')
        if a1==None or a2==None or at==None:
            missing += [['at', t]]
            c = 8*[0.0]
        else:
            if len(at) == 3: at += at
            c = at + [a1[0], a2[0]]
        coeffs.append(tuple(c))
    return coeffs,missing
