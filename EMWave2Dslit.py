# EMWave2D.py
# Simulasi gel. EM 2D
#
# using class Picture
#
# I Wayan Sudiarta
# Fisika, FMIPA, Universitas Mataram
# last updated: 14 Oktober 2018
#
# % python EMWave2D.py
#

import stddraw
import math
from picture import Picture
from color import Color

MAX_GRAY_SCALE = 255
pic = Picture(201,201)
stddraw.setCanvasSize(pic.width(), pic.height())
width  = pic.width()
height = pic.height()
i, j = 0, 0
n = 0
dx = 0.04 # dx = dy
dt = 0.3*dx # untuk stabil dt < dx/2^1/2
dtdx = dt/dx
   
ex = [[0.0] * 201 for i in range(201)]
ey = [[0.0] * 201 for i in range(201)]
hz = [[0.0] * 201 for i in range(201)]
   
# Materi
NM = 2
media = [[0] * 201 for i in range(201)]
epsR = [0.0]*NM
epsI = [0.0]*NM
muR = [0.0]*NM
muI = [0.0]*NM
s = [0.0]*NM
sm = [0.0]*NM
ch1 = [0.0]*NM
ch2 = [0.0]*NM
ce1 = [0.0]*NM
ce2 = [0.0]*NM
  
omega, freq = 0.0, 1.0
hue = 0.0
   
# nilai permitivitas relatif
epsR[0] = 1.0
epsI[0] = 0.0

epsR[1] = 2
epsI[1] = 0.0 
	  
muR[0] = 1.0
muI[0] = 0.0

muR[1] = 1.0
muI[1] = 0.0

omega = 2*math.pi*freq
	  
s[0] = omega*epsI[0]
s[1] = omega*epsI[1]
	  
ch1[0] = (1.0 - 0.5*dt*sm[0]/muR[0])/(1.0 + 0.5*dt*sm[0]/muR[0])
ch1[1] = (1.0 - 0.5*dt*sm[1]/muR[1])/(1.0 + 0.5*dt*sm[1]/muR[1])

ch2[0] = dtdx/(muR[0]*(1.0 + 0.5*dt*sm[0]/muR[0]))
ch2[1] = dtdx/(muR[1]*(1.0 + 0.5*dt*sm[1]/muR[1]))
	  
ce1[0] = (1.0 - 0.5*dt*s[0]/epsR[0])/(1.0 + 0.5*dt*s[0]/epsR[0])
ce1[1] = (1.0 - 0.5*dt*s[1]/epsR[1])/(1.0 + 0.5*dt*s[1]/epsR[1])

ce2[0] = dtdx/(epsR[0]*(1.0 + 0.5*dt*s[0]/epsR[0]))
ce2[1] = dtdx/(epsR[1]*(1.0 + 0.5*dt*s[1]/epsR[1]))
	  
# initialisasi media 0
for i in range(201):
    for j in range(201):
        media[i][j] = 0

#	medium 1  
#for(i=200 i<=400 i++){
#    for(j=0 j<=400 j++){
#      media[i][j] = 1

		
for i in range(201):
    for j in range(201):
        ex[i][j] = 0.0 
        ey[i][j] = 0.0 
        hz[i][j] = 0.0 

n = 0
while 1:
    n = n + 1
    # Hitung Hz
    for i in range(1, 200):
        for j in range(1, 200):
            hz[i][j] = ch1[media[i][j]]*hz[i][j] + ch2[media[i][j]]*(ex[i][j+1] - ex[i][j] - ey[i+1][j] + ey[i][j])

    # Hitung Ex
    for i in range(1, 200):
        for j in range(1, 200):
            ex[i][j] = ce1[media[i][j]]*ex[i][j] + ce2[media[i][j]]*(hz[i][j] - hz[i][j-1])
    
    # dinding ada celah
    for i in range(1, 200):
        if i < 60 or i > 140:
            ex[i][100] = 0.0
        if i > 80 and i < 120:
            ex[i][100] = 0.0

    # Hitung Ey
    for i in range(1, 200):
        for j in range(1, 200):
            ey[i][j] = ce1[media[i][j]]*ey[i][j] + ce2[media[i][j]]*(hz[i-1][j] - hz[i][j])

    hz[100][50] += math.sin(omega*n*dt)

    if n%20 == 0:
        stddraw.clear()
        for i in range(201):
            for j in range(201):
                v = (MAX_GRAY_SCALE / 2.0)  + (500*hz[i][j])
                if v < 0:
                    grayScale = 0
                elif v > MAX_GRAY_SCALE:
                    grayScale = MAX_GRAY_SCALE
                else:
                    grayScale = int(v)            
                color = Color(grayScale, grayScale, grayScale)
                pic.set(i, j, color)

        stddraw.picture(pic)
        stddraw.show(20)
