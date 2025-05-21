import numpy as np


def valores_propios_rossler(a, b, c):
    # Ajustar la h a partir de los valores propios
    if ((c**2)-(4*a*b) >= 0):
        z_eq = (c - np.sqrt((c**2)-(4*a*b)))/(2*a)
        x_eq = a*z_eq
        y_eq = -z_eq
    else:
        x_eq = 0
        y_eq = 0
        z_eq = 0

    J = np.array([[0, -1, -1], [1, a, 0], [z_eq, 0, x_eq - c]])
    eigenvalues, _ = np.linalg.eig(J)

    vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
    vp_non_zero = vp != 0
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    # print(vp[vp_non_zero],vp_inverse)
    # t_step=0.001
    t_step = np.round(vp_min/150, 5)
    transient = int(vp_max*5/t_step)
    steady_state = int(500*vp_max)

    return x_eq, y_eq, z_eq, t_step

# Ejemplo de uso
a, b, c = 0.458987,	0.476328,	2.51389

x_eq, y_eq, z_eq, t_step = valores_propios_rossler(a,b,c)
print(f'Ancho de paso: h = {t_step}')
print(f'Puntos de equilibrio: x_eq = {x_eq}, y_eq = {y_eq}, z_eq = {z_eq}')
