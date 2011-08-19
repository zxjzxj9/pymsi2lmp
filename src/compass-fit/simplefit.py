#!/usr/bin/python
from scipy import *
from scipy.optimize import leastsq
import matplotlib.pyplot as py

# This is what the leastsq function will try to minimize.
def residuals(p, y, x):
    return y-residuals.f(x,p)

# Pairwise interaction (see LAMMPS pair style, class II)
def vdw(x, p):
    return p[0]*(2.0*(p[1]/x)**9.0 - 3.0*(p[1]/x)**6.0)

# Quadratic bond (see LAMMPS bond style, class II)
def bond(x, p):
    dx = x - p[0]
    return p[1]*dx**2 + p[2]*dx**3 + p[3]*dx**4

# Used for fitting pairwise (van Der Wall) interactions.
if False:
    # For c3"
    x = [2.0,2.5,3.0,3.5,4.0]
    y = [41.625859,4.236509,0.430629,-0.028536,-0.063023]
    # For H1n
    x = [0.75,    1.0,      1.5,       2.0,      2.5]
    y = [5.98067, 0.287862, -0.00974, -0.00325, -0.00099]
    # For n3mh
    x = [2.5,     2.75,     3.0,      3.5,       3.75]
    y = [7.321947,2.220073, 0.529387,-0.177866, -0.19909]
    residuals.f = vdw
    p = [1.0, 3.0]

# Used for fitting simple bonded interactions.
if True:
    # For bond c3" n3mh
#    x = [1.0,        1.5,      2.0,        2.25,       1.75]
#    y = [146.983661, 4.588004, 174.828352, 582.655124, 42.892288]
    # For bond c3a n3mh
#    x = [1.0,        1.5,      2.0,        2.25,       1.75]
#    y = [118.754598,3.162439,118.471945,390.187052,30.418417]
    # For bond o1= c3"
    #y = [62.383171, 38.393136,471.256458,1436.898517,136.129059]
    # For bond n3mh h1n
#    y = [0.047344, 76.257574,916.053145,56.864239,289.965334]
#    x = [1.0,1.5,2,0.75,1.75]
    # For bond c3' n3mh
#    x = [1.0,        1.5,      2.0,        2.25,       1.75]
#    y = [110.293263, 4.101014, 169.932014, 549.07061,  41.039037]


    # For angle n3mh c3' o1=
#    x = [100.0,  110.0,  120.0,   130.0,    140.0]
#    y = [14.694994, 4.020057, 0.067717, 1.901066, 8.334059]

    #for bond c4o c4o
    x = [1.75, 1.0, 1.5,    1.25,    2.0]
    y = [13.544727, 180.126874, 0.0,29.225039, 54.684375]




    residuals.f = bond
    p = [1.5, 1.0, 0.0, 0.0]

plsq = leastsq(residuals, p, xtol=1e-10,args=(y, x), maxfev=2000)
dx = max(x) - min(x)
xx = linspace(min(x)-0.05*dx, max(x)+0.05*dx, 300)
yy = [ residuals.f(xi, plsq[0]) for xi in xx ]
for xi,yi in zip(x,y): print xi, abs((residuals.f(xi, plsq[0])-yi)/yi)

py.plot(xx,yy,'r', x, y, 'bx')
py.show()

print plsq[0]
