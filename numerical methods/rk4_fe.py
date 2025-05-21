import numpy as np
import matplotlib.pyplot as plt
# =============================================================
#Método de Runge Kutta de 4to orden
def runge_kutta4(y,h):
    k1 = function(y)
    k2 = function(y + h*0.5*np.array(k1))
    k3 = function(y + h*0.5*np.array(k2))
    k4 = function(y + h*np.array(k3))
    y = y + (h/6.)*(np.array(k1) + 2*np.array(k2) + 2*np.array(k3) + np.array(k4))
    return y
# =============================================================
# Ecuaciones del sistema de Lorenz
def function(y):
    dy = y
    return [dy]
# =============================================================
num_steps = 100
h = 0.1
t_start = 0
t_end = h*num_steps
y0 = 1

sol_Euler = np.zeros((num_steps+2, 1))
sol_RK4 = np.zeros((num_steps+2, 1))

sol_Euler[0] = [y0]
sol_RK4[0] = [y0]

# Resolver el sistema de Lorenz usando Forward Euler
for i in range(num_steps+1):
    dy = function(sol_Euler[i])
    sol_Euler[i+1] = sol_Euler[i] + h * np.array([dy])
# =============================================================
# Resolver el sistema de Lorenz usando Runge-Kutta 4th
for i in range(num_steps+1):
    sol_RK4[i+1] = runge_kutta4(sol_RK4[i],h)

for i in range(num_steps+1):
    print(h*i,sol_Euler[i],sol_RK4[i])

# Crear un vector de tiempo
t = np.linspace(t_start, t_end, num_steps+2)

# Calcular la solución analítica
sol_analitica = np.exp(t)

# Graficar las soluciones
plt.figure(figsize=(10,6))
plt.plot(t, sol_Euler, label='Euler', color='blue')
plt.plot(t, sol_RK4, label='Runge-Kutta 4', color='red', linestyle='--')
plt.plot(t, sol_analitica, label='Solución Analítica', color='green', linestyle='-.')
plt.xlabel('Tiempo')
plt.ylabel('y(t)')
plt.title('Comparación de Métodos: Euler vs Runge-Kutta 4 vs Solución Analítica')
plt.legend()
plt.grid(True)
plt.show()
