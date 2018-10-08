#  EMWave1D.py
#  Simulation EM Wave 1D using FDTD method  
#
#  I Wayan Sudiarta
#  Fisika, Fakultas MIPA, Universitas Mataram
#  last updated: 8 Oktober 2018
#
#  Maxwell's Equation 1D
#  Plane wave propagating in +x direction
#  linear polarized in Y direction 
#
#  dEy/dt = - dBz/dx
#  dBz/dt = - dEy/dx
#  
#  Using Yee's grid scheme:
#  bz^n+1/2[i+1/2] = bz^n-1/2[i+1/2] - dtdivdx*(ey^n[i+1] - ey^n[i])
#  ey^n+1[i] = ey^n[i] - dtdivdx*(bz^n+1/2[i+1/2] - bz^n+1/2[i-1/2])
# 
#  run:  % python EMWave1D.py
#

import stddraw
import math

# set the canvas size
stddraw.setCanvasSize(500,250);

# set the axis scales
stddraw.setXscale(0.0, 10.0);
stddraw.setYscale(-2.5, 2.5);

i = 0
n = 0
NX = 500
ic = 250
dx = 0.02

# For stable simulation dt < dx
dt = dx/2.0; 
dtdivdx = dt/dx

ey = [0.0]*(NX+1)
bz = [0.0]*(NX+1)

# position for plotting Ey and Bz
x1, y1, x2, y2 = 0.0, 0.0, 0.0, 0.0
freq = 1.0
 
omega = 2.0*math.pi*freq
	  
# (1) initialize array to zero
for i in range(NX+1):
    ey[i] = 0.0 
    bz[i] = 0.0

n=0
while True:
    n = n + 1
    # (2) Update Bz
    for i in range(1,NX):
        bz[i] = bz[i] - dtdivdx*(ey[i+1] - ey[i])

    # (3) Generate an EM source
    # (a) a sinusoidal wave
    # bz[ic] = bz[ic] + math.sin(omega*n*dt);
    # (b) an EM pulse
    if n<60:
        bz[ic] = bz[ic] + 2.0*math.exp(-0.005*(n-30.0)*(n-30.0))
    
    # (4) Update Ey
    for i in range(1,NX):
        ey[i] = ey[i] - dtdivdx*(bz[i] - bz[i-1])

    # (5) Draw Ey fields
    stddraw.clear()
    stddraw.setPenColor(stddraw.RED) 
    for i in range(NX):
        x1 = i*dx
        y1 = ey[i]
        x2 = (i+1)*dx
        y2 = ey[i+1]
        stddraw.line(x1, y1, x2, y2) 
    
    # Display and wait for 20 ms
    stddraw.show(20);

