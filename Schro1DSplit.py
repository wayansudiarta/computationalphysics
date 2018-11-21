#  Schro1DSplit.py
#  Simulation Time-Dependent Schrodinger 1D 
#  Using FDTD method  
#
#  I Wayan Sudiarta
#  Fisika, Fakultas MIPA, Universitas Mataram
#  last updated: 21 November 2018
#
#  Split Schrodinger Equation 1D
#  A particle propagating in +x direction
#
#  dpsiR/dt = - (1/2)d^2 psiI/dx^2 + V(x) psiI
#  dpsiI/dt =   (1/2)d^2 psiR/dx^2 - V(x) psiR
#
#  Ref:  P.B. Visscher, A Fast Explicit Algorithm 
#        for the Time-Dependent Schrodinger Equation,
#        Computers In Physics, 596–598 (Nov/Dec 1991).  
#
#  run:  % python Schro1DSplit.py
#

import stddraw
import cmath
import math

# set the canvas size
stddraw.setCanvasSize(500,250);

# set the axis scales
stddraw.setXscale(0.0, 20.0);
stddraw.setYscale(-2.5, 2.5);

i = 0
n = 0
NX = 200
ic = 100
dx = 20.0/NX

# For stable simulation dt < dx
dt = dx*dx/2.0; 
cc = 0.5*dt/(dx*dx)

psiR = [0.0]*(NX+1)
psiI = [0.0]*(NX+1)
v = [0.0]*(NX+1)

# position for plotting Ey and Bz
x1, y1, x2, y2 = 0.0, 0.0, 0.0, 0.0
freq = 1.0
 
# (1) initialize array to zero
for i in range(NX+1):
    psiR[i] = 0.0 
    psiI[i] = 0.0
    v[i] = 0.0

# initial values for psiR and psiI
p0 = 4.0
p02 = p0*p0
dp = 1
dp2 = dp*dp

for i in range(1,NX-1):
    # psiR at t = 0
    t = 0.0;
    x = (i-ic)*dx
    ta = complex(-dp2*x*x/2.0, -p02*t/2.0 + p0*x)
    tb = complex(1.0, dp2*t)
    ta = ta/tb
    ta = cmath.exp(ta)
    tc = complex(dp,0.0)
    tc = tc/tb
    tc = cmath.sqrt(tc)
    tc = tc*ta  		  
    psiR[i]= tc.real
    # psiI at times dt/2
    t = dt/2.0
    x = (i-ic)*dx
    ta = complex(-dp2*x*x/2.0, -p02*t/2.0 + p0*x)
    tb = complex(1.0, dp2*t)
    ta = ta/tb
    ta = cmath.exp(ta)
    tc = complex(dp,0.0)
    tc = tc/tb
    tc = cmath.sqrt(tc)
    tc = tc*ta  		  
    psiI[i] = tc.imag


n=0
while True:
    n = n + 1
    # (2) Update psiR
    for i in range(1,NX-1):
        psiR[i] = psiR[i] - cc*(psiI[i+1] - 2*psiI[i] + psiI[i-1]) + v[i]*psiI[i]*dt

    # (4) Update psiI
    for i in range(1,NX-1):
        psiI[i] = psiI[i] + cc*(psiR[i+1] - 2*psiR[i] + psiR[i-1]) - v[i]*psiR[i]*dt


		# (5) Draw psiR 
    stddraw.clear()
    stddraw.setPenColor(stddraw.RED) 
    for i in range(NX):
        x1 = i*dx
        y1 = psiR[i]
        x2 = (i+1)*dx
        y2 = psiR[i+1]
        stddraw.line(x1, y1, x2, y2)
				
    stddraw.setPenColor(stddraw.BLUE) 
    for i in range(NX):
        x1 = i*dx
        y1 = psiI[i]
        x2 = (i+1)*dx
        y2 = psiI[i+1]
        stddraw.line(x1, y1, x2, y2) 

    stddraw.setPenColor(stddraw.GREEN) 
    for i in range(NX):
        x1 = i*dx
        y1 = math.sqrt(psiR[i]*psiR[i] + psiI[i]*psiI[i])
        x2 = (i+1)*dx
        y2 = math.sqrt(psiR[i+1]*psiR[i+1] + psiI[i+1]*psiI[i+1])
        stddraw.line(x1, y1, x2, y2) 
    
    # Display and wait for 5 ms
    stddraw.show(5);

