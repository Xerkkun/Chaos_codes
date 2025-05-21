import sys
import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt

###############################################################################
# PARÁMETROS Y CONSTANTES EFORK
###############################################################################

# if len(sys.argv) != 2:
#     print("Uso: python lorenz_frac.py N1")
#     sys.exit(1)


#N1 = int(sys.argv[1])  # número de pasos de tiempo
N1 = 10000
T = 100.0               # tiempo total de simulación
alpha = 0.9935           # orden fraccionario
h = T / N1             # paso de integración
ha = h**alpha          # h^alpha

# Parámetro para "memoria corta"
# Ej: Lm = 5.0 (tiempo real), se tomarán solo los últimos (Lm/h) pasos de memoria
Lm = 10.0
nu = int(Lm / h)  # número de pasos que abarca la memoria

c2 = ( 1.0/( 2*gamma(1+alpha) ))**(1.0/alpha)
c3 = ( 1.0/( 4*gamma(1+alpha) ))**(1.0/alpha)

a21 = 1.0/(2*gamma(alpha+1)**2)
a31 = ( gamma(alpha+1)**2 * gamma(2*alpha+1) + 2*gamma(2*alpha+1)**2 -
        gamma(3*alpha+1) )/( 4*gamma(alpha+1)**2 *
        ( 2*gamma(2*alpha+1)**2 - gamma(3*alpha+1) ) )
a32 = - gamma(2*alpha+1)/( 4*( 2*gamma(2*alpha+1)**2 - gamma(3*alpha+1) ) ) 

w1 = ( 8.0*gamma(1+alpha)**3 * gamma(1+2*alpha)**2 
       - 6.0*gamma(1+alpha)**3 * gamma(1+3*alpha) 
       + gamma(1+2*alpha)*gamma(1+3*alpha) ) \
     / ( gamma(1+alpha)*gamma(1+2*alpha)*gamma(1+3*alpha) )

w2 = 2.0*gamma(1+alpha)**2 *( 4.0*gamma(1+2*alpha)**2 - gamma(1+3*alpha) ) / \
     ( gamma(1+2*alpha)*gamma(1+3*alpha) )

w3 = -8.0*gamma(1+alpha)**2 *( 2.0*gamma(1+2*alpha)**2 - gamma(1+3*alpha) )/ \
     ( gamma(1+2*alpha)*gamma(1+3*alpha) )

###############################################################################
# SISTEMA DE LORENZ FRACCIONARIO
###############################################################################

sigma = 10.0
rho   = 28.0
beta  = 8.0 / 3.0

def lorenz_system(x, y, z):
    dx = sigma*(y - x)
    dy = x*(rho - z) - y
    dz = x*y - beta*z
    return np.array([dx, dy, dz])

###############################################################################
# MEMORIA FRACCIONARIA
###############################################################################

vtn = np.zeros(N1+1)
vxn = np.zeros(N1+1)
vyn = np.zeros(N1+1)
vzn = np.zeros(N1+1)

# Condiciones iniciales
x0, y0, z0 = -10, -10, 30

vtn[0] = 0.0
vxn[0] = x0
vyn[0] = y0
vzn[0] = z0

def memory_fractional(k, t, arr):
    """
    Computa el término de memoria para una variable 'arr' (vxn, vyn o vzn).
    Se hace la sumatoria desde max(0, k-nu) hasta k-1, en lugar de 0..k-1,
    usando el 'short memory principle'.
    """
    # Límite inferior de la ventana de memoria
    start_idx = max(0, k - nu)
    
    v = 0.0
    for i in range(start_idx, k):
        t0 = vtn[i]
        t1 = vtn[i+1]
        v1 = (t - t0)**(1.0 - alpha)
        v2 = (t - t1)**(1.0 - alpha)
        a  = (arr[i+1] - arr[i]) / (h * gamma(2.0 - alpha))
        v += a*(v1 - v2)
    return v

def fn_lorenz(k, t, x, y, z):
    dxdt, dydt, dzdt = lorenz_system(x, y, z)
    mem_x = memory_fractional(k, t, vxn)
    mem_y = memory_fractional(k, t, vyn)
    mem_z = memory_fractional(k, t, vzn)

    fx = dxdt - mem_x
    fy = dydt - mem_y
    fz = dzdt - mem_z
    return np.array([fx, fy, fz])

###############################################################################
# PRIMER PASO (n = 0)
###############################################################################

n = 0
tn = 0.0
x_n, y_n, z_n = x0, y0, z0

# K1, K2, K3 para el EFORK en el primer paso (sin memoria aún)
K1 = ha * lorenz_system(x_n, y_n, z_n)
K2 = ha * lorenz_system(x_n + a21*K1[0],
                        y_n + a21*K1[1],
                        z_n + a21*K1[2])
K3 = ha * lorenz_system(x_n + a31*K2[0] + a32*K1[0],
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

###############################################################################
# BUCLE PARA n > 0
###############################################################################
while n < N1-1:
    n += 1
    tn = n*h
    
    # K1, K2, K3 con memoria incluida
    K1 = ha * fn_lorenz(n, tn, x_n, y_n, z_n)
    K2 = ha * fn_lorenz(n, tn + c2*h,
                        x_n + a21*K1[0],
                        y_n + a21*K1[1],
                        z_n + a21*K1[2])
    K3 = ha * fn_lorenz(n, tn + c3*h,
                        x_n + a31*K2[0] + a32*K1[0],
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


# Especifica la ruta completa donde se guardará el archivo
ruta = r'D:\INAOE\Doctorado\STM32\EFORK_short.txt'  # Cambia esta ruta según tus necesidades

# Combina las soluciones en una matriz, donde cada columna representa t, x, y, z
datos = np.column_stack((vtn, vxn, vyn, vzn))

# Guarda la matriz en el archivo de texto usando un delimitador de tabulador.
# La opción header añade una fila con nombres de columnas y comments='' evita comentarios automáticos.
np.savetxt(ruta, datos, delimiter='\t', header='t\tx\ty\tz', comments='')


print(f"Soluciones guardadas en: {ruta}")
###############################################################################
# GRÁFICAS
###############################################################################

print("Simulación finalizada. Graficando resultados...")

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(vxn, vyn, vzn)
ax.set_title(f'Sistema de Lorenz fraccionario, alpha={alpha}, Lm={Lm}')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()

plt.figure()
plt.plot(vtn, vxn, label='x(t)')
plt.plot(vtn, vyn, label='y(t)')
plt.plot(vtn, vzn, label='z(t)')
plt.title(f'Variables de Lorenz fraccionario, alpha={alpha}, Lm={Lm}')
plt.xlabel('t')
plt.ylabel('valor')
plt.grid(True)
plt.legend()
plt.show()
