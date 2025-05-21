import os
import numpy as np
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


def runge_kutta4(X, h, sigma, rho, beta):
    k1 = lorenz_equations(X, sigma, rho, beta)
    k2 = lorenz_equations(X + h*0.5*np.array(k1), sigma, rho, beta)
    k3 = lorenz_equations(X + h*0.5*np.array(k2), sigma, rho, beta)
    k4 = lorenz_equations(X + h*np.array(k3), sigma, rho, beta)
    X = X + (h/6.)*(np.array(k1) + 2*np.array(k2) +
                    2*np.array(k3) + np.array(k4))
    return X


def mod255(x, y, z, sel):
    var2 = (x, y, z)
    bin = format(int((var2[sel]*1E6) % 255), 'b')
    while len(bin) < 8:
        bin = '0' + bin
    return bin


def valores_propios_lorenz(sigma, rho, beta):
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
    t_step = np.round(vp_min/10, 5)
    tt = int(vp_max*5/t_step)
    return x_eq, y_eq, z_eq, t_step, tt


# ===============================================================================
# Inicio
sigma, rho, beta = 10, 28, 8./3.
x_eq, y_eq, z_eq, t_step, tt = valores_propios_lorenz(sigma, rho, beta)

nn, ss, v, b, nstp = '1e6 10 z mod255 1'.split()
n, s, t, nstep = int(float(nn)), int(ss), int(tt), int(nstp)

k = 8

while (n + t) % k != 0:
    t += 1

nt = (n + t) * s
x0, y0, z0 = np.zeros(s, dtype=float), np.zeros(
    s, dtype=float), np.zeros(s, dtype=float)
c = (-x_eq, y_eq, z_eq)

for j in range(s):
    x0[j] = (random.random()-0.5) * c[0]
    y0[j] = (random.random()-0.5) * c[1]
    z0[j] = (random.random()-0.5) * c[2]

var = ['x', 'y', 'z']
sel = var.index(v)

# Crear la carpeta si no existe
output_folder = 'output_folder'
os.makedirs(output_folder, exist_ok=True)

for file_index in range(10):  # Crear 10 archivos
    file_path = os.path.join(
        output_folder, f'{b}_{v}_lorenz_optimizado_{file_index}.bin')
    with open(file_path, 'wb') as arch:
        r, i = -1, -1
        x, y, z = 0, 0, 0

        while r < s:
            i += 1
            if i == 0:
                r += 1
                if r >= s:  # Evitar desbordamiento del índice
                    break
                xn, yn, zn = x0[r], y0[r], z0[r]
            else:
                xn, yn, zn = x, y, z

            x, y, z = runge_kutta4([xn, yn, zn], t_step, sigma, rho, beta)

            if x > 1E3 or x < -1E3:
                x0[r] = (random.random() - 0.5) * c[0]
                y0[r] = (random.random() - 0.5) * c[1]
                z0[r] = (random.random() - 0.5) * c[2]
                i, r = -1, r - 1

            if i > int(t/k) - 1 and (i % nstep) == 0:
                bin = mod255(x, y, z, sel)
                arch.write(bin.encode())

            if i == ((nstep * n + t) / k) - 1:
                if r < s - 1:
                    arch.write(b'\n')
                i = -1
