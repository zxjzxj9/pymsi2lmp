#!/usr/bin/python

"""
  run_msi2lmp.py - This module provides an interface to the code that converts
  .car and .mdf data from Material Studio.
"""
import sys
sys.path.append(sys.path[0] +'/../pymsi2lmp')
import pymsi2lmp
from modify_frc import term_parameter_count as count
         
# Calls pymsi2lmp and returns the which parameters are missing and 
# the total number of unknown coefficients.
def call_msi2lmp(model, frc):
    missing = pymsi2lmp.msi2lmp(model, frc)    
    unknown = {}
    for m in missing:
        if not m[0] in unknown: unknown[m[0]] = []
        unknown[m[0]].append(m[1])
    ct = sum([sum([count(m,a) for a in unknown[m]]) for m in unknown])
    return unknown,ct