import numpy as np
from numerical_methods_v2 import *
from graficas import grafica3d

def main():
    name = 'lorenz'
    method = 'GL'
    d = 3
    
    alpha = 0.992
    
    # Definición del intervalo de tiempo y el paso de integración
    t_start = 0.0
    t_end = 100
    t_step = 0.01
    h_alpha = t_step**alpha
    
    # Definición del número de puntos en la suma de Grünwald-Letnikov
    num_steps = int((t_end-t_start)/t_step) 

    # Longitud de memoria
    Lm = 300.0

    # Número de coeficientes binomiales
    m = int(Lm/t_step)

    #Principio de memoria corta
    if num_steps < m:
        nu,mm = 1,num_steps
    else:
        nu,mm = m,m
    
    # Inicialización del vector de estado y tiempo
    x = np.zeros((num_steps+1, d))
    t = np.linspace(t_start, t_end, len(x))
    x_t = np.zeros((num_steps+1,1))
    x_t[:,0] = t
    
    # Condición inicial
    x[0,:] = 0.1, 0.1, 0.1

    sigma = 10
    rho = 28
    beta = 8/3
    
    # T = 6
    # R = 20

    # m = 0
    # n = 0

    #x, status = forward_euler(x,x_t, num_steps, t_step, m, n, d, name)
    # x, status = runge_kutta4(x, x_t, num_steps, t_step, m, n, d, name)
    x, status = grunwald_letnikov(x,h_alpha,num_steps,alpha,x_t,nu,d,mm,sigma,rho,beta,name)
    #x, status = caputo_derivative(x,x_t, num_steps, t_step, m, n, d, name)
    #x, status = forward_euler(x,x_t, num_steps, t_step, T, R, d, name)
    
    # for i in range(1,num_steps+1):
    #     dy = lorenz_equations(x[i-1,:],sigma,rho,beta)
    #     x[i,:] = x[i-1,:] + t_step * dy
    
    t, xx, y, z = x_t[:,0], x[:,0], x[:,1], x[:,2]
    grafica3d(xx,y,z,t,name,method)

if __name__=='__main__':
    main()