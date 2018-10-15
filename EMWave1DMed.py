#  EMWave1DMed.py
#  Simulation EM Wave 1D di dalam Medium  
#
#  I Wayan Sudiarta
#  Fisika, Fakultas MIPA, Universitas Mataram
#  last updated: 10 Oktober 2018
#
#  Maxwell's Equation 1D
#  Plane wave propagating in +x direction
#  linear polarized in Y direction 
#
#  dBz/dt = - dEy/dx
#  dEy/dt = - 1/(mur*epsr)dBz/dx -sigma/epsr Ey
#  
#  Using Yee's grid scheme:
#  bz^n+1/2[i+1/2] = bz^n-1/2[i+1/2] - dtdivdx*(ey^n[i+1] - ey^n[i])
#  ey^n+1[i] = ce1[i] ey^n[i] - ce2[i]*(bz^n+1/2[i+1/2] - bz^n+1/2[i-1/2])
# 
#  run:  % python EMWave1DMed.py
#

import stddraw
import math

# set the canvas size
stddraw.setCanvasSize(1000,500);

# set the axis scales
stddraw.setXscale(0.0, 10.0);
stddraw.setYscale(-2.5, 2.5);

i = 0
n = 0
NX = 500
ic = 250
dx = 0.02

NM = 2  

# For stable simulation dt < dx
dt = dx/2.0; 
dtdivdx = dt/dx

ey = [0.0]*(NX+1)
bz = [0.0]*(NX+1)

epsR = [0.0]*NM
epsI = [0.0]*NM
muR = [0.0]*NM
sigma = [0.0]*NM
ce1 = [0.0]*NM
ce2 = [0.0]*NM

#indeks medium
med = [0]*(NX+1)

# position for plotting Ey and Bz
x1, y1, x2, y2 = 0.0, 0.0, 0.0, 0.0
freq = 1.0
 
omega = 2.0*math.pi*freq
  
# (1) initialize array to zero
for i in range(NX+1):
    ey[i] = 0.0 
    bz[i] = 0.0
    med[i] = 0

nr = 1.5 # ref index of water
ni = 0.0

# setting parameters
# medium 0 
epsR[0]  = 1.0
epsI[0]  = 0.0
muR[0]   = 1.0 
sigma[0] = 0.0

# medium 1
epsR[1]  = nr*nr - ni*ni
epsI[1]  = 2*nr*ni
muR[1]   = 1.0
sigma[1] = omega*epsI[1]

for i in range(NM):
    temp = 1 + 0.5*sigma[i]*dt/epsR[i]
    ce1[i] = (1 - 0.5*sigma[i]*dt/epsR[i])/temp
    ce2[i] = dt/(muR[i]*epsR[i]*dx*temp)


for i in range(NX+1):
    if i < (NX/2): 
        med[i] = 0
    elif i < (3*NX/4):
        med[i] = 1
    else:
        med[i] = 0

n=0
while True:
    n = n + 1
    # (2) Update Bz
    for i in range(0,NX-1):
        bz[i] = bz[i] - dtdivdx*(ey[i+1] - ey[i])

    # (3) Generate an EM source
    # (a) a sinusoidal wave
    # bz[ic - 100] = bz[ic - 100] + math.sin(omega*n*dt)
    # (b) an EM pulse
    if n<60:
        bz[ic - 100] = bz[ic - 100] + 2.0*math.exp(-0.005*(n-30.0)*(n-30.0))
    
    # (4) Update Ey
    for i in range(1,NX):
        ey[i] = ce1[med[i]]*ey[i] - ce2[med[i]]*(bz[i] - bz[i-1])

    # (5) Draw Ey fields
    stddraw.clear()
    stddraw.setPenColor(stddraw.YELLOW) 
    stddraw.filledRectangle(5,-2.5,2.5,5)
    stddraw.setPenColor(stddraw.BLUE)
    stddraw.text(2, 2, "Medium 0")
    stddraw.text(7, 2, "Medium 1")

    stddraw.setPenColor(stddraw.RED) 
    for i in range(NX):
        x1 = i*dx
        y1 = ey[i]
        x2 = (i+1)*dx
        y2 = ey[i+1]
        stddraw.line(x1, y1, x2, y2) 
    
    # Display and wait for 20 ms
    stddraw.show(10);

