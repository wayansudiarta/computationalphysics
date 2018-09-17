# bandul.py
# Simulasi gerak bandul 
# dengan metode Verlet versi posisi
#
# by I Wayan Sudiarta
# Physics, Universitas Mataram
# Mataram, NTB, Indonesia
#
# Rumus:
# xn+1 = 2xn - xn-1 + an dt^2
# x1 = x0 + v0 dt + 1/2 a0 dt^2
# supaya tidak mengubah kode, digunakan
# x = sudut theta
# v = kecepatan sudut
# a = percepatan sudut = -(g/l) sin(theta)
# -----------------------------------------

from math import pi
from math import sin, cos

n, x0 = 0, 0.0
v0 = 0.0
a0 = 0.0
an = 0.0
xn = 0.0
xnp1 = 0.0
xnm1 = 0.0
tn = 0.0
dt = 0.1
dt2 = dt*dt
g = 9.81
L = 1.0

# sudut awal dalam radian
x0 = 20*pi/180.0
v0 = 0.0
a0 = -(g/L)*sin(x0)

# nilai awal
xnm1 = x0
xn = x0 + v0*dt + 0.5*a0*dt2

for n in range(2, 1000):
    tn = n*dt

    # hitung an
    an = -(g/L)*sin(xn)

    # hitung xnp1
    xnp1 = 2*xn - xnm1 + an*dt2

    # simpan nilai
    xnm1 = xn
    xn = xnp1

    print(tn, xn)

