import numpy as np
from numerical_methods_v2 import *
from graficas import temporal_evolution
from graficas import graph_colors
from graficas import graph_colors_3d
from graficas import grafica3d

def main():
    name = 'rossler'
    method = 'RK4'
    dd = 3
    # Lorenz
    # a, b, c = 10, 28, 8/3
    # a, b, c = 5.301300510626442808e+00 , 1.352804249862039967e+02 , 7.840922407503264635e-01
    alpha = 0.992
    # Rossler
    a, b, c = 0.2, 0.2, 5.7

    # =============================================================
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
    # Verificar que si se genere caos
    
    vp = np.array( [eigenvalues.real,eigenvalues.imag], dtype=float )
    vp_non_zero = vp != 0 
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    t_step = np.round(vp_min/10, 5)

    transient = int(vp_max*5/t_step)
    steady_state = int(1E6) #int(500*vp_max) #int(1E5)
    num_steps = transient+steady_state
    n = int(1E4)

    print('Ancho de paso: ', t_step)
    print('Número de pasos: ', num_steps)
    print('Transitorio:, ', transient)
    print('Estacionario:, ', steady_state)
    # =============================================================
    # Definición del intervalo de tiempo y el paso de integración
     # Tiempo de integración
    t_start = 0
    t_end = num_steps*t_step    
    # ==============================================================
    # Solo Grünwald-Letnikov =======================================
    h_alpha = t_step**alpha
    
    # Definición del número de puntos en la suma de Grünwald-Letnikov
    # num_steps = int((t_end-t_start)/t_step) 

    # Longitud de memoria
    Lm = 300.0

    # Número de coeficientes binomiales
    m = int(Lm/t_step)

    #Principio de memoria corta
    if num_steps < m:
        nu,mm = 1,num_steps
    else:
        nu,mm = m,m
    # =============================================================
    # Inicialización del vector de estado y tiempo
    x = np.zeros((num_steps+1, dd))
    t = np.linspace(t_start, t_end, len(x))
    x_t = np.zeros((num_steps+1,1))
    x_t[:,0] = t
    
    # Condición inicial
    x[0, :] = 0.12434987362700446, 0.6717225928431056, 0.7588531905022086

    # x[0,:] = x_eq, -y_eq, z_eq
    print('Condiciones iniciales: ', x[0,:])    

    #Wei-Chen
    # m, n = 10, 100, -0.1966

    #Barati 
    # m, n = 0.6, 0.3, 0.5
    
    
    
    #Moore-Spiegel
    # T, R = 6, 20

    #Wang
   # m, n = 0, 0
    
    # Rossler
    # m, n = 0.2, 0.2, 5.7
    
    #chen
    # m, n = 35, 3, 28

    #x, status = forward_euler(x,x_t, num_steps, t_step, m, n, d, name)
    x, status = runge_kutta4(x, x_t, num_steps, t_step, a, b, c, dd, name)
    #x, status = grunwald_letnikov(x,h_alpha,num_steps,alpha,x_t,nu,d,mm,m,n,name)
    #x, status = caputo_derivative(x,x_t, num_steps, t_step, m, n, d, name)
    #x, status = forward_euler(x,x_t, num_steps, t_step, T, R, d, name)
    
    t, xx, y, z = x_t[:,0], x[:,0], x[:,1], x[:,2]
    
    temporal_evolution(xx, y, z, t, t_step, transient*t_step,name, method)
    graph_colors(xx, y, z, t, name, method)
    graph_colors_3d(xx, y, z, t, name, method)
    grafica3d(xx,y,z,t,name,method)
if __name__=='__main__':
    main()
