# leapfrog.py
# Simulasi gerak jatuh bebas 
# dengan metode leapfrog
# 
# by I Wayan Sudiarta
# Physics, Universitas Mataram
# Mataram, NTB, Indonesia
#
# Rumus:
# xn+1 = xn + vn+1/2 dt
# vn+1/2 = vn-1/2 + an dt
# ------------------------------

n, x0 = 0, 0.0
v0 = 5.0
a0, an = -9.81, -9.81
xn = 0.0
vn = 0.0
t0, tn = 0.0, 0.0
dt = 0.1

# hitung v1/2
vn = v0 + 0.5*a0*dt
xn = x0

for n in range(1, 1000):
    tn = t0 + n*dt

    # xn+1
    xn = xn + vn*dt

    # hitung an+1
    an = -9.81

    # vn+1/2
    vn = vn + an*dt

    print(tn, vn, xn)

    if xn < 0:
        break
