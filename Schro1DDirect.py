#  Schro1DDirect.py
#  Simulation Time-Dependent Schrodinger 1D 
#  Using FDTD method  
#
#  I Wayan Sudiarta
#  Fisika, Fakultas MIPA, Universitas Mataram
#  last updated: 21 November 2018
#
#  A particle propagating in 1D, x axis 
#  i d psi/dt = - (1/2)d^2 psi/dx^2 + V(x) psi
#  Discretization: See ref. Askar and Cakmak 
#  J. Chem. Phys. 68, 2794 (1978); 
#  https://doi.org/10.1063/1.436072
#  
#	 run:  % python Schro1DDirect.py
#

import stddraw
import cmath

# set the canvas size
stddraw.setCanvasSize(500,250);

# set the axis scales
stddraw.setXscale(0.0, 20.0);
stddraw.setYscale(-0.5, 1.5);

i = 0
n = 0
NX = 400
ic = 200
dx = 20.0/NX

# For stable simulation dt < dx
dt = dx*dx/2.0; 
cc = dt/(dx*dx)

psi = [complex(0,0)]*(NX+1)   # psin   
psip1 = [complex(0,0)]*(NX+1) # psinp1
psim1 = [complex(0,0)]*(NX+1) # psinm1 
v = [0.0]*(NX+1)

# position for plotting psi
x1, y1, x2, y2 = 0.0, 0.0, 0.0, 0.0
 
# (1) initialize array to zero
for i in range(NX+1):
    psi[i] = 0.0*1j 
    psip1[i] = 0.0*1j
    psim1[i] = 0.0*1j
    v[i] = 0.0

# initial values for psiR and psiI
p0 = 4.0
p02 = p0*p0
dp = 0.5
dp2 = dp*dp

for i in range(1,NX-1):
    # psiR at t = 0
    t = -dt;
    x = (i-ic)*dx
    ta = complex(-dp2*x*x/2.0, -p02*t/2.0 + p0*x)
    tb = complex(1.0, dp2*t)
    ta = ta/tb
    ta = cmath.exp(ta)
    tc = complex(dp,0.0)
    tc = tc/tb
    tc = cmath.sqrt(tc)
    tc = tc*ta  		  
    psim1[i]= tc
    # psi at time dt
    t = 0
    x = (i-ic)*dx
    ta = complex(-dp2*x*x/2.0, -p02*t/2.0 + p0*x)
    tb = complex(1.0, dp2*t)
    ta = ta/tb
    ta = cmath.exp(ta)
    tc = complex(dp,0.0)
    tc = tc/tb
    tc = cmath.sqrt(tc)
    tc = tc*ta  		  
    psi[i] = tc


n=0
while True:
    n = n + 1
    # (2) Update psi
    for i in range(1,NX-1):
        psip1[i] = psim1[i] + cc*1j*(psi[i+1] - 2*psi[i] + psi[i-1]) - 2*1j*v[i]*psi[i]*dt

    # (4) Save
    for i in range(1,NX-1):
        psim1[i] = psi[i]
        psi[i] = psip1[i]
				
		# (5) Draw psiR 
    stddraw.clear()
    stddraw.setPenColor(stddraw.RED) 
    for i in range(NX):
        x1 = i*dx
        y1 = abs(psi[i])
        x2 = (i+1)*dx
        y2 = abs(psi[i+1])
        stddraw.line(x1, y1, x2, y2)
				
    # Display and wait for 5 ms
    stddraw.show(5);

