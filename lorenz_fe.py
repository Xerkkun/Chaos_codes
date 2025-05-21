import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def lorenz_equations(X, sigma, rho, beta):
    x, y, z = X
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]

def forward_euler(x0, t_start, t_end, t_step, sigma, rho, beta):
    num_steps = int((t_end - t_start) / t_step)
    x = np.zeros((num_steps+1, 3))
    x[0] = x0

    for i in range(num_steps):
        dy = lorenz_equations(x[i], sigma, rho, beta)
        x[i+1] = x[i] + t_step * np.array(dy)

    return x

def main():
    # Parámetros del sistema de Lorenz
    sigma = 10.0
    rho = 28.0
    beta = 8.0/3.0

    # Condiciones iniciales
    x0 = [1.0, 1.0, 1.0]

    # Definición del intervalo de tiempo y el paso de integración
    t_start = 0.0
    t_end = 50.0
    t_step = 0.01

    x = forward_euler(x0, t_start, t_end, t_step, sigma, rho, beta)

    # Graficar
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x[:,0], x[:,1], x[:,2])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Lorenz Attractor using Forward Euler')
    plt.show()

if __name__ == '__main__':
    main()
