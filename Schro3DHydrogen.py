#  Schro3DHydrogen.py
#  Simulation Time-Dependent Schrodinger 3D 
#  Using FDTD method  
#
#  I Wayan Sudiarta
#  Fisika, Fakultas MIPA, Universitas Mataram
#  last updated: 19 December 2018
#
#  Schrodinger Equation 3D in imaginary time
#  An electron in atom hydrogen 
#
#  dpsi/dt = (1/2)d^2 psi/dx^2 - V(r) psi
#
#  V(r) = -1/r
#
#  For FDTD Computation potential is shifted to
# 
#  V(r) = -1/r + 1/dx
#
#  Ref:  
#  I. W. Sudiarta, D. J. W. Geldart
#  Journal of Physics A: Mathematical and Theoretical 40 (8), 1885
#
#
#  run:  % python Schro3DHydrogen.py
#

import math
# use NumPy for array
import numpy as np
import stddraw

# set the canvas size
stddraw.setCanvasSize(500,500);

# set the axis scales
stddraw.setXscale(-3, 3);
stddraw.setYscale(0, 1);

i = 0
n = 0

NX = 40

ic = NX//2
jc = NX//2
kc = NX//2
 
dx = 6.0/NX
dx2 = dx*dx
dx3 = dx2*dx

coef = 0;

print(dx)
print(1/dx)

# For stable simulation dt < dx/3
dt = 0.5*(dx2/3.0); 

psi = np.zeros((NX+1, NX+1, NX+1))
psiNew = np.zeros((NX+1, NX+1, NX+1))
a = np.zeros((NX+1, NX+1, NX+1))
b = np.zeros((NX+1, NX+1, NX+1))
v = np.zeros((NX+1, NX+1, NX+1))

# position for plotting Ey and Bz
x1, y1, x2, y2 = 0.0, 0.0, 0.0, 0.0
 
# (1) initialize array to zero
for i in range(NX+1):
   for j in range(NX+1):
       for k in range(NX+1):
           # zeros
           psi[i, j, k] = 0.0 
           psiNew[i, j, k] = 0.0

           # potential 			  
           x = (i - ic)*dx
           y = (j - jc)*dx
           z = (k - kc)*dx
           r2 = x*x + y*y + z*z
           r = math.sqrt(r2)
           if r < dx:
               v[i, j, k] = 0.0
           else:
               v[i, j, k] = -1.0/r + 1.0/dx
  
# initial psi  and coeffs    
for i in range(1, NX-1):
   for j in range(1, NX-1):
       for k in range(1, NX-1):
           psi[i, j, k] = 1.0 
           temp = 1.0 + 0.5*dt*v[i, j, k]    
           a[i, j, k] = (1.0 - 0.5*dt*v[i, j, k])/temp
           b[i, j, k] = dt/(2.0*dx2*temp)

n=0
while True:
    n = n + 1
    # (2) Update psi
    for i in range(1, NX-1):
        for j in range(1, NX-1):
            for k in range(1, NX-1):
                temp  = psi[i+1, j, k] + psi[i-1, j, k]
                temp += psi[i, j+1, k] + psi[i, j-1, k]
                temp += psi[i, j, k+1] + psi[i, j, k-1]
                temp -= 6.0*psi[i, j, k]
                psiNew[i, j, k] = a[i, j, k]*psi[i, j, k] + b[i, j, k]*temp

    # (3) Save
    for i in range(1, NX-1):
        for j in range(1, NX-1):
            for k in range(1, NX-1):
                psi[i, j, k] = psiNew[i, j, k]

    # (3) Hitung Energi setiap 10 iterasi
    if (n%10) == 0:
        integral1 = 0.0
        integral2 = 0.0
        for i in range(1, NX-1):
            for j in range(1, NX-1):
                for k in range(1, NX-1):
                    integral1 += psi[i, j, k]*psi[i, j, k]
                    d2psidx2  = -0.5*(psi[i+1, j, k] - 2*psi[i, j, k] + psi[i-1, j, k])/dx2
                    d2psidx2 += -0.5*(psi[i, j+1, k] - 2*psi[i, j, k] + psi[i, j-1, k])/dx2
                    d2psidx2 += -0.5*(psi[i, j, k+1] - 2*psi[i, j, k] + psi[i, j, k-1])/dx2
                    integral2 += psi[i, j, k]*(d2psidx2 + v[i, j, k]*psi[i, j, k])
        energi = integral2/integral1 - 1.0/dx
        print(energi)
        stddraw.clear()
        stddraw.setPenColor(stddraw.BLUE) 
        for i in range(NX):
            x1 = (i-ic)*dx
            y1 = psi[i, jc, kc]/psi[ic, jc, kc]
            x2 = (i+1 - ic)*dx
            y2 = psi[i+1, jc, kc]/psi[ic,jc,kc]
            stddraw.line(x1, y1, x2, y2)
        stddraw.show(5)



