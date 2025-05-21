import numpy as np


def valores_propios_chen(a, b, c):
    # Ajustar la h a partir de los valores propios
    x_eq = np.sqrt(2*b*c-a*b)
    y_eq = x_eq
    z_eq = np.power(x_eq, 2)/b

    J = np.array([[-a, a, 0], [c-a-z_eq, c, -x_eq], [y_eq, x_eq, -b]])
    eigenvalues, _ = np.linalg.eig(J)

    vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
    vp_non_zero = vp != 0
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    # print(vp[vp_non_zero],vp_inverse)
    # t_step=0.001
    t_step = np.round(vp_min/50, 5)
    transient = int(vp_max*5/t_step)
    steady_state = int(500*vp_max)

    return x_eq, y_eq, z_eq, t_step


# Ejemplo de uso
a, b, c = 2.821888891274968714e+01,  7.066636808327453334e-01,  2.475217691457946856e+01
x_eq, y_eq, z_eq, t_step = valores_propios_chen(a, b, c)
print(f'Ancho de paso: h = {t_step}')
print(f'Puntos de equilibrio: x_eq = {x_eq}, y_eq = {y_eq}, z_eq = {z_eq}')
