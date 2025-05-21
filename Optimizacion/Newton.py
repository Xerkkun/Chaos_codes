import numpy as np
import math
import sys

n = len(sys.argv)
if n!=2:
  print("Args: valor_inicial_x valor_inicial_y")
  sys.exit(1)

x = float(sys.argv[1])
y = float(sys.argv[2])

#print(vx.shape)

J = np.zeros( (2,1) )
H = np.zeros( (2,2) )
invH = np.zeros( (2,2) )

a = 1
b = 100

def Jacobiana(J,x,a,b):
    J[0,0] = 2*(x-a)-4*b*x(y-x*x)
    J[1,0] = 2*b*(y-x*x)

def Hessiana(H,x,a,b):
    H[0,0] = 2-4*b*y+12*b*x*x
    H[0,1] = -4*b*x
    H[1,0] = -4*b*x
    H[1,1] = 2*b

i=0
while i<20 :
    Jacobiana(J,x,a,b)
    Hessiana(H,x,a,b)

    det = H[0,0]*H[1,1] - H[0,1]*H[1,0]
    invH[0,0] = H[1,1]/det
    invH[0,1] = -H[0,1]/det
    invH[1,0] = -H[1,0]/det
    invH[1,1] = H[0,0]/det

    #Calcular la inversa de H
    vdx = -invH @ J

    n = vdx[0,0]*vdx[0,0] + vdx[1,0]*vdx[1,0]
    n = math.sqrt( n )

    if n < 1e-10:
        break

    x += vdx[0,0]
    y += vdx[1,0]

    print( i,x,y)

    i += 1
