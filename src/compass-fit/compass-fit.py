#!/usr/bin/env python
"""

    SYNOPSIS
   
      compass-fit.py [basename] [original_frc_path]
"""

from run_msi2lmp import * 
import sys
import modify_frc
from scipy.optimize import leastsq
import lammps
from numpy import *
from numpy.linalg import norm
from discover_output import read_arc_file, read_disco_energy_file

# Main input file.
def main(args):
    import glob
    # Finds the first *.mdf file in the directory and uses that as the basename.
    base    = glob.glob('*.mdf')[0][0:-4]
    frcfile = '../../data/compass.frc'
    # If a base name is specified then override the automatically found one.
    if len(args) > 1: base     = args[1]
    if len(args) > 2: frcfile  = args[2]
    if len(args) > 3: 
        print 'compass_fitter.py [basename] [original_frc_path]'
        return

    # Calls msi2lmp and returns the output as strings.
    run_lmp.missing, param_count = call_msi2lmp(base, frcfile)
    run_lmp.base    = base    
    # Initialize unknown parameters
    v = zeros((param_count))
    print 'Attempting to fit',param_count,'parameters.'        
    
    residual.x, residual.e = read_arc_file(base)        
    if param_count > 0:                
        # Minimize error in energy with least squares.
        v = leastsq(residual, v, args=(), epsfcn=0.0004)[0]        
            
    err = norm(residual(v)) / norm(residual.e)     
    print 'Relative error norm is: %.2g' % err 
    
    # Last run energy breakdown.
    lmp_e = lammps.get_pe_from_log_file()
    ms_e  = read_disco_energy_file(base)
    
    def compare(lmp, ms):
        if ms != 0.0:
            return '%11.7f\t%11.7f\t%11.7f\t%5.3g' %(lmp, ms, lmp-ms,100.0*abs((lmp-ms)/ms))
        else:
            return '%11.7f\t%11.7f\t%11.7f\t%5.3g' %(lmp, ms, lmp-ms,0.0)
    print '       LAMMPS             Mat Studio      Error       % difference'
    print 'pe  ', compare(lmp_e['PotEng'], ms_e['pe'])
    print 'vdw ', compare(lmp_e['E_vdwl'], ms_e['vdw'])
    print 'bond', compare(lmp_e['E_bond'], ms_e['b'])
    print 'angl', compare(lmp_e['E_angle'], ms_e['a']+ms_e['ba']+ms_e['bb'])
    print 'tor ', compare(lmp_e['E_dihed'], ms_e['tor']
                                           +ms_e['mbt']
                                           +ms_e['ebt']
                                           +ms_e['at']
                                           +ms_e['aat']
                                           +ms_e['bb13'])
    print 'oop ', compare(lmp_e['E_impro'], ms_e['oop']+ms_e['aa'])
    print 'coul', compare(lmp_e['E_coul'], ms_e['q'])   
                 
# Returns the difference between the MS and LAMMPS energies.
def residual(v):
    energy = run_lmp(v, residual.x)
    print norm(energy - residual.e) / norm(residual.e) 
    return energy - residual.e

    
# Returns the LAMMPS energies for a set of trajectories.
def run_lmp(v, x):
    energy = zeros((len(x)))
    # Makes a new frc file in the current directory with the new parameters.
    modify_frc.update('../../data/compass.frc', 'compass1.frc', run_lmp.missing, v)
    call_msi2lmp(run_lmp.base, 'compass1.frc')
    for i, xi in enumerate(x):
        lammps.modify_data_file(run_lmp.base, xi, i) 
        lmpout    = lammps.run(run_lmp.base, i, '/opt/lammps/lmp_openmpi')
        energy[i] = lammps.extract_energy(lmpout)
    return array(energy)

if __name__ == '__main__': main(sys.argv)

