import numpy as np


def valores_propios_lorenz(sigma, rho, beta):
    # Calcular los puntos de equilibrio
    x_eq, y_eq, z_eq = np.sqrt(
        np.abs(beta*(rho-1))), np.sqrt(np.abs(beta*(rho-1))), rho-1

    # Calcular la matriz Jacobiana en el punto de equilibrio
    J = np.array(
        [[-sigma, sigma, 0],
         [rho - z_eq, -1, -x_eq],
         [y_eq, x_eq, -beta]]
    )

    # Calcular los valores propios de la matriz Jacobiana
    eigenvalues, _ = np.linalg.eig(J)

    # Preparar los resultados
    vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
    vp_non_zero = vp != 0
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    t_step = np.round(vp_min/10, 5)
    tt = int(vp_max*5/t_step)

    return x_eq, y_eq, z_eq, t_step


# Ejemplo de uso
sigma, rho, beta = 5.301300510626442808e+00, 1.352804249862039967e+02, 7.840922407503264635e-01

x_eq, y_eq, z_eq, t_step = valores_propios_lorenz(sigma, rho, beta)
print(f'Ancho de paso: h = {t_step}')
print(f'Puntos de equilibrio: x_eq = {x_eq}, y_eq = {y_eq}, z_eq = {z_eq}')
