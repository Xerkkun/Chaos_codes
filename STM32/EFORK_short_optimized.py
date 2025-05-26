import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt

###############################################################################
# PARÁMETROS Y CONSTANTES EFORK
###############################################################################

N1 = 50000              # número de pasos de tiempo
T = 500.0               # tiempo total de simulación
alpha = 0.985         # orden fraccionario
h = T / N1              # paso de integración
ha = h**alpha           # h^alpha

# Parámetro para "memoria corta"
Lm = 1.0
nu = int(Lm / h)        # número de pasos que abarca la memoria

c2 = (1.0/(2*gamma(1+alpha)))**(1.0/alpha)
c3 = (1.0/(4*gamma(1+alpha)))**(1.0/alpha)

a21 = 1.0/(2*gamma(alpha+1)**2)
a31 = ( gamma(alpha+1)**2 * gamma(2*alpha+1) + 2*gamma(2*alpha+1)**2 - gamma(3*alpha+1) ) / ( 4*gamma(alpha+1)**2 * (2*gamma(2*alpha+1)**2 - gamma(3*alpha+1)) )
a32 = - gamma(2*alpha+1) / ( 4*(2*gamma(2*alpha+1)**2 - gamma(3*alpha+1)) )

w1 = (8.0*gamma(1+alpha)**3 * gamma(1+2*alpha)**2 - 6.0*gamma(1+alpha)**3 * gamma(1+3*alpha) + gamma(1+2*alpha)*gamma(1+3*alpha)) / (gamma(1+alpha)*gamma(1+2*alpha)*gamma(1+3*alpha))
w2 = 2.0*gamma(1+alpha)**2 * (4.0*gamma(1+2*alpha)**2 - gamma(1+3*alpha)) / (gamma(1+2*alpha)*gamma(1+3*alpha))
w3 = -8.0*gamma(1+alpha)**2 * (2.0*gamma(1+2*alpha)**2 - gamma(1+3*alpha)) / (gamma(1+2*alpha)*gamma(1+3*alpha))

###############################################################################
# SISTEMA DE LORENZ FRACCIONARIO
###############################################################################
# Parametros del sistema de Rossler
a = 0.2
b = 0.2 
c = 5.7

# Parametros del sistema de Lorenz
sigma = 10.0
rho   = 28.0
beta  = 8.0 / 3.0

# Parametros el sistema de Chen
u = 7.5
v = 1.0
w = 5.0

def lorenz_system(x, y, z):
    """Calcula el lado derecho del sistema de Lorenz"""
    dx = sigma*(y - x)
    dy = x*(rho - z) - y
    dz = x*y - beta*z
    return np.array([dx, dy, dz])

def chen_system(x, y, z):
    """Calcula el lado derecho del sistema de Chen"""
    dx = u*(y - x)
    dy = (w - u)*x - x*z + w*y
    dz = x*y - v*z
    return np.array([dx, dy, dz])

def rossler_system(x, y, z):
    """Calcula el lado derecho del sistema de Rossler"""
    dx = -y - z
    dy = x + a*y
    dz = b + z*(x - c)
    return np.array([dx, dy, dz])

###############################################################################
# FUNCIONES PARA EL MÉTODO DE GRÜNWALD-LETNIKOV
###############################################################################
# Versión vectorizada del término de memoria
def memory_fractional(k, t, arr, vtn, h, alpha, nu):
    """
    Computa el término de memoria para una variable 'arr' (vxn, vyn o vzn)
    de forma vectorizada, eliminando el bucle for interno.
    """
    start_idx = max(0, k - nu)
    idx = np.arange(start_idx, k)
    t0 = vtn[idx]
    t1 = vtn[idx+1]
    v1 = (t - t0)**(1.0 - alpha)
    v2 = (t - t1)**(1.0 - alpha)
    denom = h * gamma(2.0 - alpha)
    a = (arr[idx+1] - arr[idx]) / denom
    return np.sum(a * (v1 - v2))

def fn_rossler(k, t, x, y, z, vtn, vxn, vyn, vzn, h, alpha, nu):
    """
    Calcula la función del sistema de rossler incluyendo el término de memoria.
    """
    dxdt, dydt, dzdt = rossler_system(x, y, z)
    mem_x = memory_fractional(k, t, vxn, vtn, h, alpha, nu)
    mem_y = memory_fractional(k, t, vyn, vtn, h, alpha, nu)
    mem_z = memory_fractional(k, t, vzn, vtn, h, alpha, nu)
    fx = dxdt - mem_x
    fy = dydt - mem_y
    fz = dzdt - mem_z
    return np.array([fx, fy, fz])

###############################################################################
# SIMULACIÓN DEL SISTEMA
###############################################################################

# Inicialización de las variables (memoria para t, x, y, z)
vtn = np.zeros(N1+1)
vxn = np.zeros(N1+1)
vyn = np.zeros(N1+1)
vzn = np.zeros(N1+1)

# Condiciones iniciales
x0, y0, z0 = 1, 0, 0 
vtn[0] = 0.0
vxn[0] = x0
vyn[0] = y0
vzn[0] = z0

# Primer paso (n = 0) sin incluir memoria
n = 0
tn = 0.0
x_n, y_n, z_n = x0, y0, z0

K1 = ha * rossler_system(x_n, y_n, z_n)
K2 = ha * rossler_system(x_n + a21*K1[0],
                        y_n + a21*K1[1],
                        z_n + a21*K1[2])
K3 = ha * rossler_system(x_n + a31*K2[0] + a32*K1[0],
                        y_n + a31*K2[1] + a32*K1[1],
                        z_n + a31*K2[2] + a32*K1[2])

x_n1 = x_n + w1*K1[0] + w2*K2[0] + w3*K3[0]
y_n1 = y_n + w1*K1[1] + w2*K2[1] + w3*K3[1]
z_n1 = z_n + w1*K1[2] + w2*K2[2] + w3*K3[2]

tn1 = (n+1)*h
vtn[n+1] = tn1
vxn[n+1] = x_n1
vyn[n+1] = y_n1
vzn[n+1] = z_n1

x_n, y_n, z_n = x_n1, y_n1, z_n1

# Bucle principal para n > 0
for n in range(1, N1):
    tn = n * h

    # Se calculan K1, K2 y K3 incluyendo el término de memoria
    K1 = ha * fn_rossler(n, tn, x_n, y_n, z_n, vtn, vxn, vyn, vzn, h, alpha, nu)
    K2 = ha * fn_rossler(n, tn + c2*h,
                         x_n + a21*K1[0],
                         y_n + a21*K1[1],
                         z_n + a21*K1[2],
                         vtn, vxn, vyn, vzn, h, alpha, nu)
    K3 = ha * fn_rossler(n, tn + c3*h,
                         x_n + a31*K2[0] + a32*K1[0],
                         y_n + a31*K2[1] + a32*K1[1],
                         z_n + a31*K2[2] + a32*K1[2],
                         vtn, vxn, vyn, vzn, h, alpha, nu)

    x_n1 = x_n + w1*K1[0] + w2*K2[0] + w3*K3[0]
    #y_n1 = y_n + w1*K1[1] + w2*K2[1] + w3*K2[1] if False else None  # <- Nota: corregir error tipográfico si ocurre
    y_n1 = y_n + w1*K1[1] + w2*K2[1] + w3*K3[1]
    z_n1 = z_n + w1*K1[2] + w2*K2[2] + w3*K3[2]

    tn1 = (n+1)*h
    vtn[n+1] = tn1
    vxn[n+1] = x_n1
    vyn[n+1] = y_n1
    vzn[n+1] = z_n1

    x_n, y_n, z_n = x_n1, y_n1, z_n1

# Guardar los resultados en un archivo de texto
ruta = r'C:\Users\moren\Desktop\Chaos_codes\Chaos_codes\STM32\EFORK_rossler_python.txt'
datos = np.column_stack((vtn, vxn, vyn, vzn))
np.savetxt(ruta, datos, delimiter='\t', header='t\tx\ty\tz', comments='')
print(f"Soluciones guardadas en: {ruta}")

###############################################################################
# GRÁFICAS
###############################################################################
print("Simulación finalizada. Graficando resultados...")

"""Graficar atractores en 3D, guardar en PDF y mostrar la gráfica."""
plt.figure (figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.plot(vxn, vyn, "m", lw=0.3)
plt.xlabel("x")
plt.ylabel("y")

plt.subplot(1, 3, 2)
plt.plot(vxn, vzn, "m", lw=0.3)
plt.xlabel("x")
plt.ylabel("z")

plt.subplot(1, 3, 3)
plt.plot(vyn, vzn, "m", lw=0.3)
plt.xlabel("y")
plt.ylabel("z")

plt.tight_layout()
plt.savefig("rossler_atractores_EFORK_python.pdf", dpi=400, bbox_inches='tight')
plt.show()  # Se muestra la gráfica al final
plt.clf()


# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.plot(vxn, vyn, vzn)
# ax.set_title(f'Sistema de Lorenz fraccionario, alpha={alpha}, Lm={Lm}')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# plt.show()

# plt.figure()
# plt.plot(vtn, vxn, label='x(t)')
# plt.plot(vtn, vyn, label='y(t)')
# plt.plot(vtn, vzn, label='z(t)')
# plt.title(f'Variables de Lorenz fraccionario, alpha={alpha}, Lm={Lm}')
# plt.xlabel('t')
# plt.ylabel('valor')
# plt.grid(True)
# plt.legend()
# plt.show()
