#  atomLJ2.py
#  Simulasi 2 atom LJ pada ruang 2D
#  
#  I Wayan Sudiarta
#  Fisika, FMIPA, Universitas Mataram
#
#  Dependencies: stddraw.py
#
#
#  % python atomLJ2.py
#
#

import stddraw
import math

# Setting ukuran layar
stddraw.setXscale(-3, 3)
stddraw.setYscale(-3, 3)

# variabel indeks
i, j, k = 0, 0, 0

# Deklarasi Array posisi kecepatan percepatan
# untuk dua atom LJ
# Array [0.0, 0.0]
x = [0.0]*2
y = [0.0]*2
vx = [0.0]*2
vy = [0.0]*2
ax = [0.0]*2
ay = [0.0]*2

dx, dy = 0.0, 0.0
r, ri, ri3, ri6 = 0.0, 0.0, 0.0, 0.0
force = 0.0	

# Variabel Waktu
n = 0
dt = 0.01
dt2 = dt*dt

# Jari-jari lingkaran
radius = 0.4

# Energi kinetik dan potensial
ke, pe = 0.0, 0.0

# Nilai awal posisi, kecepatan dan percepatan
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

# loop
n = 0
while True: 
    n = n + 1
    # update menggunakan metode Verlet versi kecepatan

    # step(1) update posisi
    for i in range(2):
        x[i] = x[i] + vx[i]*dt + 0.5*ax[i]*dt2
        y[i] = y[i] + vy[i]*dt + 0.5*ay[i]*dt2 

    # step(2) update kecepatan 1/2
    for i in range(2):
        vx[i] = vx[i] + 0.5*ax[i]*dt
        vy[i] = vy[i] + 0.5*ay[i]*dt

    # step(3) Update percepatan
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    r = math.sqrt(dx*dx + dy*dy)
    ri = 1/r
    ri3 = ri*ri*ri
    ri6 = ri3*ri3
    force = 24.0*ri*ri6*(2.0*ri6 - 1.0)/r

    ax[0] = -force*dx
    ay[0] = -force*dy
    ax[1] = force*dx
    ay[1] = force*dy
  
     
    # step(4) update kecepatan 1/2
    ke = 0;
    for i in range(2):
        vx[i] = vx[i] + 0.5*ax[i]*dt
        vy[i] = vy[i] + 0.5*ay[i]*dt

        # compute kinetic energy
        ke = ke + 0.5*(vx[i]*vx[i] + vy[i]*vy[i])

    # hitung energi kinetik dan potensial
    ke = 0;
    for i in range(2):
        ke = ke + 0.5*(vx[i]*vx[i] + vy[i]*vy[i])

    pe = 4*ri6*(ri6 - 1.0);
    print(ke, "  ", pe, "  ", (pe+ke))

    # step(5) Visualisasi Gambar lingkaran
    stddraw.setPenColor(stddraw.BLUE)
    for i in range(2):
        stddraw.filledCircle(x[i], y[i], radius) 
    
    # display and pause for 40 ms
    stddraw.show(40); 

    # hapus lingkaran dengan warna background (WHITE)
    stddraw.setPenColor(stddraw.WHITE); 
    for i in range(2):
        stddraw.filledCircle(x[i], y[i], radius+1) 

