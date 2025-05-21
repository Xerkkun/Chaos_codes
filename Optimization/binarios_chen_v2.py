import numpy as np
import math
import random
import pandas as pd

# ===============================================================================
# Funciones
# ===============================================================================

def chen(X, aa, bb, cc):
    x, y, z = X
    dx = aa*(y - x)
    dy = (cc - aa)*x - x*z + cc*y
    dz = x*y - bb*z
    return [dx, dy, dz]

def runge_kutta4(X, h, aa, bb, cc):
    k1 = chen(X, aa, bb, cc)
    k2 = chen(X + h*0.5*np.array(k1), aa, bb, cc)
    k3 = chen(X + h*0.5*np.array(k2), aa, bb, cc)
    k4 = chen(X + h*np.array(k3), aa, bb, cc)
    X = X + (h/6.)*(np.array(k1) + 2*np.array(k2) +
                    2*np.array(k3) + np.array(k4))
    return X


def umbral(x, y, z, sel):
    u = (0, 0, 0, 0)
    var2 = (x, y, z)
    if var2[sel] > u[sel]:
        bin = '1'
    else:
        bin = '0'
    return bin


def mod255(x, y, z, sel):
    var2 = (x, y, z)
    bin = format(int((var2[sel]*1E6) % 255), 'b')
    while len(bin) < 8:  # Acompleta los numeros binarios a palabras de 8 bits
        bin = '0' + bin
    return bin


def valores_propios_chen(aa, bb, cc):
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
    t_step = np.round(vp_min/10, 5)

    tt = int(vp_max*5/t_step)
    return x_eq, y_eq, z_eq, t_step, tt


def run_simulation(aa, bb, cc, output_file):
    x_eq, y_eq, z_eq, t_step, tt = valores_propios_chen(aa, bb, cc)

    nn, ss, v, b, nstp = input(
        'Parámetros de entrada (n, s, v, b, nstep): ').split()
    n, s, t, nstep = int(float(nn)), int(ss), int(tt), int(nstp)

    if b == "umbral":
        k = 1
    elif b == "mod255":
        k = 8
    else:
        print("Método no definido")

    while ((n+t) % k) != 0:
        t = t + 1

    nt = (n+t)*s  # pasos totales

    x0, y0, z0 = np.zeros(s, dtype=float), np.zeros(
        s, dtype=float), np.zeros(s, dtype=float)

    c = (4.246, 4.728, 13.470)

    for j in range(0, s):
        x0[j] = (random.random()-0.5)*c[0]
        y0[j] = (random.random()-0.5)*c[1]
        z0[j] = (random.random()-0.5)*c[2]

    print("Número de pasos:", n)
    print("Número de corridas:", s)
    print("Estado transitorio:", t)
    print("Variable para sec. binarias: ", v)
    print("Metodo sec. binarias: ", b)
    print('Ancho de paso: ', t_step)
    print('Parámetros de diseño: ', aa, bb, cc)

    arch = open(output_file, "wb")
    r, i = -1, -1

    var = ['x', 'y', 'z']
    sel = var.index(v)

    x, y, z = 4.246, 4.728, 13.470  # condiciones iniciales

    while r < s-1:
        i = i + 1
        if i == 0:
            r = r + 1
            xn, yn, zn = x0[r], y0[r], z0[r]
        else:
            xn, yn, zn = x, y, z

        x, y, z = runge_kutta4([xn, yn, zn], t_step, aa, bb, cc)

        if (x > 1E3 or x < -1E3):
            print("Overflow in r = ", r)
            print(x)
            antes = arch.tell()
            x0[r] = (random.random()-0.5)*c[0]
            y0[r] = (random.random()-0.5)*c[1]
            z0[r] = (random.random()-0.5)*c[2]
            print(x0[r], y0[r], z0[r])
            pos = ((n*r)+r)-antes
            arch.seek(pos, 1)
            i, r = -1, r-1

        if i > int(t/k)-1 and (i % nstep) == 0:
            if b == "umbral":
                bin = umbral(x, y, z, sel)
            elif b == "mod255":
                bin = mod255(x, y, z, sel)
            arch.write(bin.encode())

        if i == ((nstep*n+t)/k)-1:
            if (r < s-1):
                arch.write(("\n").encode())
            i = -1
    arch.close()


def main():
    # Leer los valores de sigma, rho y beta desde el archivo top_3_gwo.txt
    params_df = pd.read_csv(
        '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/top_3_pso_chen.txt', sep=" ", header=None)
    # Tomar solo las tres primeras columnas
    params_list = params_df[[0, 1, 2]].values.tolist()

    # Ejecutar la simulación para cada conjunto de parámetros
    for idx, params in enumerate(params_list):
        output_file = f'chen_bin_{idx+1}.rnd'
        run_simulation(params[0], params[1], params[2], output_file)


if __name__ == '__main__':
    main()
