import numpy as np
import math
import random
# ===============================================================================
# Funciones
# ===============================================================================


def lorenz_equations(X, sigma, rho, beta):
    x, y, z = X
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]
# ===============================================================================


def rossler(X, aa, bb, cc):
    x, y, z = X
    dx = -y - z
    dy = x + aa*y
    dz = bb + z*(x - cc)
    return [dx, dy, dz]
# ===============================================================================


def chen(X, aa, bb, cc):
    x, y, z = X
    dx = aa*(y - x)
    dy = (cc - aa)*x - x*z + cc*y
    dz = x*y - bb*z
    return [dx, dy, dz]
# ==============================================================================


def runge_kutta4(X, h, aa, bb, cc):
    k1 = rossler(X, aa, bb, cc)
    k2 = rossler(X + h*0.5*np.array(k1), aa, bb, cc)
    k3 = rossler(X + h*0.5*np.array(k2), aa, bb, cc)
    k4 = rossler(X + h*np.array(k3), aa, bb, cc)
    X = X + (h/6.)*(np.array(k1) + 2*np.array(k2) +
                    2*np.array(k3) + np.array(k4))

    return X
# ===============================================================================
# Métodos de secuencias binarias


def umbral(x, y, z, sel):
    u = (0, 0, 0, 0)
    var2 = (x, y, z)

    if var2[sel] > u[sel]:
        bin = '1'
    else:
        bin = '0'
    return bin
# -------------------------------------------------------------------------------


def mod255(x, y, z, sel):
    u = (0, 0, 0)
    var2 = (x, y, z)
    bin = format(int((var2[sel]*1E6) % 255), 'b')
    while len(bin) < 8:  # Acompleta los numeros binarios a palabras de 8 bits
        bin = '0' + bin
    return bin
# -------------------------------------------------------------------------------


def randomMap(x, y, z, sel):
    u = (0, 0, 0)
    var2 = (x, y, z)
    r, exp = math.frexp(var2[sel])

    ival = int(2251799813685248 * r)  # 2^51
    s = 0x00000000000000ff & ival
    s = s << 8

    r, exp = math.frexp(var2[sel])
    ival = int(r * 2251799813685248)
    s = s | (0x00000000000000ff & ival)
    bin = format(s, 'b')

    while len(bin) < 16:
        bin = '0' + bin
    return bin
# -------------------------------------------------------------------------------


def valores_propios_lorenz(sigma, rho, beta):

    # Ajustar la h a partir de los valores propios
    x_eq, y_eq, z_eq = np.sqrt(
        np.abs(beta*(rho-1))), np.sqrt(np.abs(beta*(rho-1))), rho-1

    J = np.array(
        [[-sigma, sigma, 0], [rho - z_eq, -1, -x_eq], [y_eq, x_eq, -beta]])
    eigenvalues, _ = np.linalg.eig(J)

    vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
    vp_non_zero = vp != 0
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    # print(vp[vp_non_zero],vp_inverse)
    # t_step=0.001
    t_step = np.round(vp_min/10, 5)

    tt = int(vp_max*5/t_step)
    return x_eq, y_eq, z_eq, t_step, tt
# -------------------------------------------------------------------------------


def valores_propios_rossler(aa, bb, cc):
    # Ajustar la h a partir de los valores propios
    if ((cc**2)-(4*aa*bb) >= 0):
        z_eq = (cc - np.sqrt((cc**2)-(4*aa*bb)))/(2*aa)
        x_eq = aa*z_eq
        y_eq = -z_eq
    else:
        x_eq = 0
        y_eq = 0
        z_eq = 0

    J = np.array([[0, -1, -1], [1, aa, 0], [z_eq, 0, x_eq - cc]])
    eigenvalues, _ = np.linalg.eig(J)

    vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
    vp_non_zero = vp != 0
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    # print(vp[vp_non_zero],vp_inverse)
    # t_step=0.001
    t_step = np.round(vp_min/150, 5)

    tt = int(vp_max*5/t_step)
    return x_eq, y_eq, z_eq, t_step, tt
# -------------------------------------------------------------------------------


def valores_propios_chen(aa, bb, cc):
    # Ajustar la h a partir de los valores propios
    x_eq = np.sqrt(2*bb*cc-aa*bb)
    y_eq = x_eq
    z_eq = np.power(x_eq, 2)/bb

    J = np.array([[-aa, aa, 0], [cc-aa-z_eq, cc, -x_eq], [y_eq, x_eq, -bb]])
    eigenvalues, _ = np.linalg.eig(J)

    vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
    vp_non_zero = vp != 0
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    # print(vp[vp_non_zero],vp_inverse)
    # t_step=0.001
    t_step = np.round(vp_min/50, 5)

    tt = int(vp_max*10/t_step)
    return x_eq, y_eq, z_eq, t_step, tt
# ===============================================================================
# Inicio

# parámetros de entrada
# n: pasos por corridas (1e6)
# s: corridas (1000)
# t: transitorio (5000)
# met: método numérico de resolución (FE,BE,RK4,AB6,AM4,G4)
# bin: método de generación de secuencias binarias (umbral,mod255)
# nstep: frecuencia de muestreo


# originales lorenz 10, 28, 8./3.
# originales rossler 0.2, 0.2, 5.7
# originales chen 35, 3, 28

# sigma, rho, beta = 10, 28, 8./3.
# sigma, rho, beta = 5.301300510626442808, 135.2804249862039967, 0.7840922407503264635
# aa, bb, cc = 0.2, 0.2, 5.7
aa, bb, cc = 0.458987,	0.476328,	2.51389


x_eq, y_eq, z_eq, t_step, tt = valores_propios_rossler(aa, bb, cc)

nn, ss, v1, v2, b, nstp = input('Parámetros de entrada: ').split()
n, s, t, nstep = int(float(nn)), int(ss), int(tt), int(nstp)

if b == "umbral":
    k = 1
elif b == "mod255":
    k = 8
elif b == "randomMap":
    k = 16
else:
    print("Método no definido")

print(t)
while ((n+t) % k) != 0:
    t = t + 1

nt = (n+t)*s  # pasos totales
print(n+t)

x0, y0, z0 = np.zeros(s, dtype=float), np.zeros(
    s, dtype=float), np.zeros(s, dtype=float)

# c = (4.246, 4.728, 13.470)
c = (1, 1, 1)
# c = (-x_eq, y_eq, z_eq)

for j in range(0, s):
    x0[j] = (random.random()-0.01)*c[0]
    y0[j] = (random.random()-0.01)*c[1]
    z0[j] = (random.random()-0.01)*c[2]

print("Número de pasos:", n)
print("Número de corridas:", s)
print("Estado transitorio:", t)
print("Variable para sec. binarias: ", v1, v2)
print("Metodo sec. binarias: ", b)
print('Ancho de paso: ', t_step)
print('Parámetros de diseño: ', aa, bb, cc)

# "wb" para escribir archivos con formato binario
arch = open(b + "_" + "_rossler_optimizado.rnd", "wb")
r, i = -1, -1

var = ['x', 'y', 'z']
sel1 = var.index(v1)
sel2 = var.index(v2)

while r < s:
    i = i + 1
    # print(i)
    if i == 0:
        r = r + 1
        xn, yn, zn = x0[r], y0[r], z0[r]
        # print(r)
    else:
        xn, yn, zn = x, y, z

    x, y, z = runge_kutta4([xn, yn, zn], t_step, aa, bb, cc)
    # print(x, y, z)

    # Condición de paro por desbordamiento
    if (x > 1E3 or x < -1E3):
        print("Overflow in r = ", r)
        print(x)
        antes = arch.tell()
        x0[r] = (random.random()-0.01)*c[0]
        y0[r] = (random.random()-0.01)*c[1]
        z0[r] = (random.random()-0.01)*c[2]
        print(x0[r], y0[r], z0[r])
        pos = ((n*r)+r)-antes
        arch.seek(pos, 1)
        i, r = -1, r-1

        # Con frecuencia de muestreo
    if i > int(t/k)-1 and (i % nstep) == 0:
        if b == "umbral":
            bin = umbral(x, y, z, sel1)

        elif b == "mod255":
            bin = mod255(x, y, z, sel1)

        elif b == "randomMap":
            bin = randomMap(x, y, z, sel1)

        arch.write(bin.encode())

    if i == ((nstep*n+t)/k)-1:
        if (r < s-1):
            arch.write(("\n").encode())
        i = -1
arch.close()
