#  molecule2.py
#  Simulasi 2 atom LJ pada ruang 2D
#  Dependencies: stddraw.py
#
#  % python molecule2.py
#

import stddraw
import math

# Setting ukuran koordinat pada layar
stddraw.setXscale(-3, 3)
stddraw.setYscale(-3, 3)

# variabel indeks
i, j, k = 0, 0, 0
m = 0

# Deklarasi Array posisi kecepatan percepatan
# array [0.0, 0.0]
x = [0.0]*2
y = [0.0]*2
vx = [0.0]*2
vy = [0.0]*2
ax = [0.0]*2
ay = [0.0]*2

dx, dy = 0.0, 0.0
r, ri, ri3, ri6 = 0.0, 0.0, 0.0, 0.0
force = 0.0	

# Waktu
n = 0
dt = 0.01
dt2 = dt*dt

radius = 0.4
ke, pe = 0.0, 0.0

# Nilai awal
m = 0
tke = 0

# partikel 0
x[0] = -1.0
y[0] = 0.0
vx[0] = 0.0
vy[0] = 0.0
ax[0] = 0.0
ay[0] = 0.0

# partikel 1
x[1] = 2.0
y[1] = 0.5
vx[1] = -2.0
vy[1] = 0.0
ax[1] = 0.0
ay[1] = 0.0

n = 0
# main loop
while True: 
    n = n + 1
    # update using Verlet method

    # update posisi
    for i in range(2):
        x[i] = x[i] + vx[i]*dt + 0.5*ax[i]*dt2
        y[i] = y[i] + vy[i]*dt + 0.5*ay[i]*dt2 

    # update kecepatan 1/2
    for i in range(2):
        vx[i] = vx[i] + 0.5*ax[i]*dt
        vy[i] = vy[i] + 0.5*ay[i]*dt

    # Update percepatan
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    # compute distance				   
    r = math.sqrt(dx*dx + dy*dy)
    ri = 1/r
    ri3 = ri*ri*ri
    ri6 = ri3*ri3
    force = 24.0*ri*ri6*(2.0*ri6 - 1.0)/r

    ax[0] = -force*dx
    ay[0] = -force*dy
    ax[1] = force*dx
    ay[1] = force*dy
  
    # compute potential energy
    pe = 4*ri6*(ri6 - 1.0);
     
    # update kecepatan 1/2
    ke = 0;
    for i in range(2):
        vx[i] = vx[i] + 0.5*ax[i]*dt
        vy[i] = vy[i] + 0.5*ay[i]*dt

        # compute kinetic energy
        ke = ke + 0.5*(vx[i]*vx[i] + vy[i]*vy[i])

    #stddraw.clear();
			
    # draw ball on the screen
    stddraw.setPenColor(stddraw.BLUE)
    for i in range(2):
        stddraw.filledCircle(x[i], y[i], radius) 
    
    # display and pause for 40 ms
    stddraw.show(40); 

    print(ke, "  ", pe, "  ", (pe+ke))

    # remove ball on the screen
    stddraw.setPenColor(stddraw.WHITE); 
    for i in range(2):
        stddraw.filledCircle(x[i], y[i], radius+1) 

