#  Schro1DDifusi.py
#  Simulation Time-Dependent Schrodinger 1D 
#  Using FDTD method  
#
#  I Wayan Sudiarta
#  Fisika, Fakultas MIPA, Universitas Mataram
#  last updated: 21 November 2018
#
#  Schrodinger Equation 1D in imaginary time
#  A particle propagating in +x direction
#
#  dpsi/dt = (1/2)d^2 psi/dx^2 - V(x) psi
#
#  Ref:  
#  I. W. Sudiarta, D. J. W. Geldart
#  Journal of Physics A: Mathematical and Theoretical 40 (8), 1885

#
#  run:  % python Schro1DDifusi.py
#

import stddraw
import cmath
import math

# set the canvas size
stddraw.setCanvasSize(500,500);

# set the axis scales
stddraw.setXscale(0.0, 1.0);
stddraw.setYscale(-1.5, 1.5);

i = 0
n = 0
NX = 100
ic = 50
dx = 1.0/NX
dx2 = dx*dx
coef = 0;

# For stable simulation dt < dx
dt = 0.5*(dx2/2.0); 

psi = [0.0]*(NX+1)
psiNew = [0.0]*(NX+1)
a = [0.0]*(NX+1)
b = [0.0]*(NX+1)
v = [0.0]*(NX+1)

# position for plotting Ey and Bz
x1, y1, x2, y2 = 0.0, 0.0, 0.0, 0.0
freq = 1.0
 
# (1) initialize array to zero
for i in range(NX+1):
    psi[i] = 0.0 
    psiNew[i] = 0.0
    v[i] = 0.0

# psi awal ?    
for i in range(1, NX-1):
    psi[i] = 1.0 
    temp = 1 + 0.5*dt*v[i]    
    a[i] = (1 - 0.5*dt*v[i])/temp
    b[i] = dt/(2*dx2*temp)

n=0
while True:
    n = n + 1
    # (2) Update psi
    for i in range(1,NX-1):
        psiNew[i] = a[i]*psi[i] + b[i]*(psi[i+1] - 2*psi[i] + psi[i-1])

    # (3) Save
    # (3) Hitung Koefisien Normalisasi
    coef = 0 
    for i in range(1,NX-1):
        coef += psi[i]*psi[i]
    coef = coef*dx
    coef = math.sqrt(1.0/coef)    

    for i in range(1,NX-1):
        psi[i] = coef*psiNew[i]
        
    if n%20 == 0:		
		# (5) Draw psi 
        stddraw.clear()
        stddraw.setPenColor(stddraw.RED) 
        for i in range(NX):
            x1 = i*dx
            y1 = psi[i]
            x2 = (i+1)*dx
            y2 = psi[i+1]
            stddraw.line(x1, y1, x2, y2)
    
        # Display and wait for 5 ms
        stddraw.show(5);

