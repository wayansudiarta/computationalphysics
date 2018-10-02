#  atomLJN.py
#  Simulasi N atom LJ pada ruang 2D
#  
#  I Wayan Sudiarta
#  Fisika, FMIPA, Universitas Mataram
#
#  Dependencies: stddraw.py
#
#  % python atomLJ2N.py
#

import stddraw
import math
import random

# Panjang box
L = 10.0
N = 36

# array untuk warna
warna = [stddraw.RED, stddraw.BLUE, stddraw.GREEN, stddraw.MAGENTA, stddraw.CYAN, stddraw.ORANGE, stddraw.PINK, stddraw.YELLOW, stddraw.GRAY, stddraw.BLACK]

# ukuran koordinat pada layar
stddraw.setXscale(0, L)
stddraw.setYscale(0, L)

# indeks
i, j, k = 0, 0, 0

# Deklarasi Array posisi kecepatan percepatan
x = [0.0]*N
y = [0.0]*N
vx = [0.0]*N
vy = [0.0]*N
ax = [0.0]*N
ay = [0.0]*N

xnew, ynew = 0.0, 0.0
dx, dy = 0.0, 0.0
r, ri, ri3, ri6 = 0.0, 0.0, 0.0, 0.0	
force = 0.0

# Waktu
n = 0
dt = 0.01
dt2 = dt*dt

# radius lingkaran
radius = L/30.0

# energi kinetik dan potensial
ke, pe = 0.0, 0.0

# Nilai awal
tke = 0

k = 0
for i in range(6):
    for j in range(6):
        x[k] = (i+1)*L/7.0
        y[k] = (j+1)*L/7.0
        k = k + 1

for i in range(N):
    vx[i] = 5*(2*random.random()-1.0)
    vy[i] = 5*(2*random.random()-1.0)
    ax[i] = 0.0;
    ay[i] = 0.0;


# loop
n = 0
while True:
    n = n + 1

    # update using Verlet method velocity

    # update posisi
    for i in range(N):
        xnew = x[i] + vx[i]*dt + 0.5*ax[i]*dt2
        ynew = y[i] + vy[i]*dt + 0.5*ay[i]*dt2 

        # Periodic Boundary Conditions
        if xnew <= 0.0:
            xnew = xnew + L
        if xnew > L:
            xnew = xnew - L
        if ynew <= 0.0:
            ynew = ynew + L
        if ynew > L:
            ynew = ynew - L

        x[i] = xnew
        y[i] = ynew

    # update kecepatan 1/2
    for i in range(N):
        vx[i] = vx[i] + 0.5*ax[i]*dt
        vy[i] = vy[i] + 0.5*ay[i]*dt

    # Update percepatan
    # set initialy to zero
    for i in range(N):
        ax[i] = 0.0
        ay[i] = 0.0

    pe = 0.0 
    for i in range(N):
        for j in range((i+1), N):
            dx = x[i] - x[j]
            dy = y[i] - y[j]

            # minimum image convention
            if math.fabs(dx) > 0.5*L:
                dx = dx - math.copysign(1, dx)*L
            if math.fabs(dy) > 0.5*L:
                dy = dy - math.copysign(1, dy)*L
 
            # compute distance				   
            r = math.sqrt(dx*dx + dy*dy)
            ri = 1/r
            ri3 = ri*ri*ri
            ri6 = ri3*ri3
            force = 24.0*ri*ri6*(2.0*ri6 - 1.0)/r

            ax[i] = ax[i] + force*dx
            ay[i] = ay[i] + force*dy
            ax[j] = ax[j] - force*dx
            ay[j] = ay[j] - force*dy
  
            # compute potential energy
            pe = pe + 4*ri6*(ri6 - 1.0)

    # update kecepatan 1/2
    ke = 0.0
    for i in range(N):
        vx[i] = vx[i] + 0.5*ax[i]*dt
        vy[i] = vy[i] + 0.5*ay[i]*dt

        # compute kinetic energy
        ke = ke + 0.5*(vx[i]*vx[i] + vy[i]*vy[i])

    # stddraw.clear();
    # draw ball on the screen
    for i in range(N):
        # pilih warna
        j = i%10
        stddraw.setPenColor(warna[j])
        stddraw.filledCircle(x[i], y[i], radius) 
    
    # display and pause for 10 ms
    stddraw.show(10); 

    print(ke, "  ", pe, "  ", (pe+ke))

    # remove ball on the screen
    stddraw.setPenColor(stddraw.WHITE); 
    for i in range(N):
        stddraw.filledCircle(x[i], y[i], radius+1) 

