import os
import math
import sys
import subprocess
import numpy as np
from pathlib import Path
from matplotlib import style
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from pymoo.optimize import minimize
from pymoo.core.problem import Problem
from scipy.fftpack import fft, fftfreq
from mpl_toolkits.mplot3d import axes3d
from pymoo.termination import get_termination
from pymoo.algorithms.soo.nonconvex.de import DE
# ======================================================================================
def binomial_coef(alpha,mm,decimal):
    """ Cálculo de los coeficientes binomiales """
    c = np.zeros((mm+1,1))
    c[0] = 1
    arch = open(name + "_coeficientes_entero" + ".rnd","w")
    for j in range(1,mm+1):
        c[j] = c[j-1] - (c[j-1]*(1+alpha)/j)

        arch.write('%.15f' % c[j] + '\n')
    arch.close()
    return c
# ======================================================================================
def grunwald_letnikov(x,h_alpha,k,alpha,x_t,nu,d,mm,decimal,sigma,rho,beta):
    """Definición del método de Grünwald-Letnikov para la derivada de orden fraccionario"""
    # Iteraciones del método
    sum_x = np.zeros(d)
    c = binomial_coef(alpha,mm,decimal)
    arch_1 = open(name + "_variables_entero" + ".rnd","w") #"wb" para escribir archivos con formato binario
    arch_2 = open(name + "_sumatoria_entero" + ".rnd","w")
    new_x = np.zeros(d)
    aux_x = np.zeros((nu, d))
    aux_x[0,:] = x[0,:]
    #===========================================================================
    for i in range(1,k+1):
        # Se calculan las sumas de los coeficientes binomiales en cada iteracion
        if nu == 1:
            for j in range(1,i+1):
                sum_x += c[j]*x[i-j,:]
            x[i,:] = ho4_equations(x[i-1,:], sigma, rho, beta)*h_alpha - sum_x
    #===========================================================================
        else:
            for j in range(1,nu+1):
                sum_x += c[j]*aux_x[nu-j,:]
            if i < mm: # Se llenan los valores de x hasta la iteracion m (maximo de memoria)
                x[i,:] = np.array(ho4_equations(x[i-1,:], sigma, rho, beta))*h_alpha - sum_x
                aux_x[i]=x[i,:]                
    #===========================================================================
           # Luego se van desechando los primeros valores de x que se habian calculado con anterioridad para ir guardando los nuevos
            # Si el  numero de iteracion (i) ya excede el valor de la memoria (m)
            # Desplazar todos los valores del vecto x una posicion hacia arriba y dejar vacio el último elemento del vector
            else:
                new_x = np.array(ho4_equations(x[i-1,:],sigma,rho,beta))*h_alpha - sum_x
                aux_x = np.concatenate((aux_x[1:],[new_x]),axis=0)
                x[i,:] = new_x        
                
        # Condición de paro por desbordamiento
        # if (x[i,0]>1E3 or x[i,0]<-1E3):
        #     status = 0
        #     print("saludos")
        #     break
        # else:
        #     status = 1
        status = 1
    #===========================================================================
        if d==3:
            #if i%100 == 0 :
                #print(x_t[i,0],x[i,0],x[i,1],x[i,2])
            arch_1.write('%.3f' % x_t[i,0] + '\t' + '%.8f' % x[i,0] + '\t' + '%.8f' % x[i,1] + '\t' + '%.8f' % x[i,2] + '\n')
            arch_2.write('%.8f' % sum_x[0] + '\t' + '%.8f' % sum_x[1] + '\t' + '%.8f' % sum_x[2] + '\n')

        elif d==4:
            #if i%1000 == 0 :
                #print(x_t[i,0],x[i,0],x[i,1],x[i,2],x[i,3])
            arch_1.write('%.3f' % x_t[i,0] + '\t' + '%.15f' % x[i,0] + '\t' + '%.15f' % x[i,1] + '\t' + '%.15f' % x[i,2] + '%.15f' % x[i,3] + '\n')
            arch_2.write('%.15f' % sum_x[0] + '\t' + '%.15f' % sum_x[1] + '\t' + '%.15f' % sum_x[2] + '%.15f' % sum_x[3] + '\n')

        sum_x = np.zeros(d)
    arch_1.close()
    arch_2.close()
    return x,status
# =============================================================
#Método de Runge-Kutta de 4to orden
def runge_kutta4(x, h, k, sigma, rho, beta):
    """ Método de integración numérica Runge-Kutta de 4to orden """
    for i in range(1,k+1):
        k1 = ho4_equations(x[i-1,:], sigma, rho, beta)
        k2 = ho4_equations(x[i-1,:] + h*0.5*np.array(k1), sigma, rho, beta)
        k3 = ho4_equations(x[i-1,:] + h*0.5*np.array(k2), sigma, rho, beta)
        k4 = ho4_equations(x[i-1,:] + h*np.array(k3), sigma, rho, beta)
        x[i,:] = x[i-1,:] + (h/6.)*(np.array(k1) + 2*np.array(k2) + 2*np.array(k3) + np.array(k4))
        
        # Condición de paro por desbordamiento
        if (x[i,0]>1E3 or x[i,0]<-1E3):
            status = 0
            #print("saludos")
            break
        else:
            status = 1
        #print(status)
    return x,status
# ======================================================================================
# Ecuaciones del sistema de Lorenz
def lorenz_equations(X, sigma, rho, beta):
    """ Ecuaciones del sistema caótico de Lorenz de orden entero """
    x, y, z = X
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]
# ======================================================================================
# Sistema con atractores ocultos
# a = sigma
# b = rho
# c = beta
def ho4_equations(X,sigma,rho,beta):
    """ Ecuaciones del sistema caótico con atractores ocultos """
    x, y, z = X
    dx = sigma*(y-x)
    dy = rho*x*z
    dz = beta*(1-x*y)
    return [dx, dy, dz]
# ======================================================================================
def lorenz_objective(X):
    """ Definir la función objetivo a optimizar (Dimension de Kaplan-Yorke) """
    
    name = "ho4"
    
    # Parámetros del sistema de Lorenz
    sigma = X[0]  #0,60
    rho = X[1]    #0.001,180
    beta = X[2]   #0.001,30
    
     # Definición del orden fraccionario
    #alpha = 0.99
    #alpha = X[3] #0.9,0.99
# ======================================================================================
    # Ajustar la h a partir de los valores propios
    x_eq, y_eq, z_eq = np.sqrt(np.abs(beta*(rho-1))),np.sqrt(np.abs(beta*(rho-1))),rho-1
    
    J = np.array([ [-sigma, sigma, 0], [rho - z_eq, -1, -x_eq], [y_eq, x_eq, -beta]   ])
    eigenvalues, _ = np.linalg.eig(J)

    vp = np.array( [eigenvalues.real,eigenvalues.imag], dtype=float )
    vp_non_zero = vp != 0 
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    t_step=np.round(vp_min/10,5)
    
    transient=int(vp_max*5/t_step)
    steady_state = int(2E4)
    num_steps=transient+steady_state
    
    # Ancho de paso fijo
    #t_step = 0.01
    transient = int(1E4)
    steady_state = int(2E4)
    num_steps = int(transient + steady_state)
          
    # print('Ancho de paso: ', t_step)
    # print('Número de pasos: ', num_steps)
    # print('Transitorio:, ', transient)
# ======================================================================================
    #Condiciones iniciales (cercanas a los puntos de equilibrio)
    #x0, y0, z0 = -x_eq, y_eq, z_eq
    #sigma = 10
    #rho = 28
    #beta = 8/3
    
    #Condiciones iniciales constantes:
    x0, y0, z0 = 0.1, 0., 0.1
# ======================================================================================
    # Tiempo de integración
    t_start = 0
    t_end = num_steps*t_step
    #num_steps = int((t_end - t_start) / t_step)
    
    # Definición del número de puntos en la suma de Grünwald-Letnikov
    k = int((t_end-t_start)/t_step) 
    
    #dimensiones del sistema
    d = 3
    
    # Inicialización del vector de estado y tiempo
    x = np.zeros((k+1, d))
    t = np.linspace(t_start, t_end, len(x))
    x_t = np.zeros((k+1,1))
    x_t[:,0] = t

#    h_alpha = t_step**alpha
    
    sol = np.zeros((num_steps+1, 3))
    sol[0] = [x0, y0, z0]
# ======================================================================================
    #Integración numérica del sistema de Lorenz con Grünwald-Letnikov
    # Definición de la parte decimal
    decimal = 16

    # Longitud de memoria
    Lm = 0.1

    # Número de coeficientes binomiales
    m = int(Lm/t_step)

    #Principio de memoria corta
    if k < m:
        nu,mm = 1,k
    else:
        nu,mm = m,m

    #sol,status = grunwald_letnikov(sol,h_alpha,k,alpha,x_t,nu,d,mm,decimal,sigma,rho,beta)
    sol,status = runge_kutta4 (sol, t_step, k, sigma, rho, beta)
# ======================================================================================
    if status == 1:
        
        # Transformada de Fourier
        
        # Quitar el transitorio
        sol2 = sol[transient:,:]
        # t2 = np.linspace(t_step*transient, t_step*num_steps, num_steps+1-transient)
        # dt2 = t2[1] - t2[0]

        Y2 = fft(sol2[:,0]) / (num_steps+1-transient)  # Transformada normalizada
        # Y3 = fft(sol2[:,1]) / (num_steps+1-transient)  # Transformada normalizada
        # Y4 = fft(sol2[:,2]) / (num_steps+1-transient)  # Transformada normalizada
        # frq2 = fftfreq(num_steps+1-transient, dt2) 
        
        sumax = sum(abs(Y2))    
        # sumay = sum(abs(Y3))    
        # sumaz = sum(abs(Y4))    
        
        # Encuentra los picos
        peaksx, _ = find_peaks(np.abs(Y2[0:len(sol2[:,0])//2]))
        # peaksy, _ = find_peaks(np.abs(Y3[0:len(sol2[:,1])//2]))
        # peaksz, _ = find_peaks(np.abs(Y4[0:len(sol2[:,2])//2]))
        num_peaksx = len(peaksx)
        # num_peaksy = len(peaksy)
        # num_peaksz = len(peaksz)
        # print(num_peaksx)
        # Encuentra picos a partir de un umbral
        # Calcular la magnitud de la transformada de Fourier
        # magnitudex = np.abs(Y2)
        # magnitudey = np.abs(Y3)
        # magnitudez = np.abs(Y4)

        # Encontrar los índices de los picos que sobrepasan un valor específico
        # thresholdx = 6
        # thresholdy = 6
        # thresholdz = 6
        # peaks_thrx, _ = find_peaks(magnitudex, height=thresholdx)
        # peaks_thry, _ = find_peaks(magnitudey, height=thresholdy)
        # peaks_thrz, _ = find_peaks(magnitudez, height=thresholdz)
        
        # num_peaks_thrx = len(peaks_thrx)
        # num_peaks_thry = len(peaks_thry)
        # num_peaks_thrz = len(peaks_thrz)
        
    #Calcular el valor más alto que alcanzan los picos
        # max_peak_value_x = np.max(magnitudex[peaksx])
        # max_peak_value_y = np.max(magnitudey[peaksy])
        # max_peak_value_z = np.max(magnitudey[peaksy])

        # print("El valor más alto que alcanzan los picos para x es", max_peak_value_x)
        # print("El valor más alto que alcanzan los picos para y es", max_peak_value_y)
        # print("El valor más alto que alcanzan los picos para z es", max_peak_value_z)

# ======================================================================================
    # Condición para descartar evoluciones temporales periódicas       
        # Valores propios conjugados con parte real positiva
        # Potencia de la transformada de Fourier mayor a un determinado valor 
        #  and max_peak_value_x > 1
        
        if(sumax > 100 and num_peaksx > 200):
            
            # Se guarda la evolucion temporal
            output_file = name + "_entero_de.dat"
            np.savetxt(output_file, sol, delimiter="  ")

            #print("Datos guardados en", output_file)

            # TISEAN
            # Se llama a Tisean para calcular el espectro de exponentes de lyapunov y la dimension KY
            # Realizar un promedio
            # parameters = (20, 30 , 50, 70, 90)
            parameters = (30, 50, 75, 150 , 250)
            DKY_sum = 0
            LE_sum = 0
            LE = 0
            
            
            for i in parameters:      
                lyap_specPATH = 'lyap_spec ' + name +'_entero_de.dat -x' + str(transient) + ' -c1,2,3 -m3,1 -k' +str(i) + ' -osalida_entero_de_' + name + '.lyaps'
                #check_output(str(lyap_specPATH), shell=True)

                # Ejecuta el comando de lyap_spec y redirige la salida estándar y la salida de error a /dev/null
                result = subprocess.run(lyap_specPATH, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                
                filename = "/media/programas/INAOE/Doctorado/Optimización/EDIESCA 2023/salida_entero_de_" + name + ".lyaps"
                if os.path.getsize(filename) == 0:
                    DKY = 1.0
                    print("vacío")
                else:
                    # Se obtiene la D-KY del archivo de salida de TISEAN
                    with open(filename, "r") as file:
                        file.seek(0,2)
                        file.seek(file.tell()-9)
                        DKY = float(file.read())
                        file.seek(0,2)
                        file.seek(file.tell()-309)
                        LE = float(file.read(12))
                        
                # Condición por si TISEAN se apendeja
                if (DKY == 3.0):
                    DKY = 1.0
                    
                DKY_sum = DKY_sum - DKY
                LE_sum = LE_sum + LE
                
            DKY_mean = DKY_sum/len(parameters)
            LE_mean = LE_sum/len(parameters)
            file2 = open("iteraciones_entero_de_" + name + ".dat", "a")
            file2.write('%.5f' % X[0] + '\t' + '%.5f' % X[1] + '\t' + '%.5f' % X[2] + '\t' + '%.5f' % DKY_mean + '\t' + '%.5f' % LE_mean + '\n')   
            file2.close()
            
            #max_peak_value_z = np.max(magnitudez[peaksz])
            #or num_peaks_thry > 5 or num_peaks_thrz > 5
            #if(max_peak_value_z > 15):
            #    DKY_mean = -1.0
                
        else: 
            DKY_mean = -1.0 
                    
        # print("Variables",X)
        # print("Potencia", sumax)
        # print("Valores propios", eigenvalues)
        # print("Dimension Kaplan-Yorke",-DKY_mean)
        
        #Graficar 
        # label_potx = 'pot=' + str(np.round(sumax,2)) + ',peaks=' + str(num_peaksx) + ',' + str(num_peaks_thrx)
        # label_poty = 'pot=' + str(np.round(sumay,2)) + ',peaks=' + str(num_peaksy) + ',' + str(num_peaks_thry)
        # label_potz = 'pot=' + str(np.round(sumaz,2)) + ',peaks=' + str(num_peaksz) + ',' + str(num_peaks_thrz)
        # label_eig = str(np.round(eigenvalues.real, 2) + np.round(eigenvalues.imag, 2) * 1j)
        # label_dky = str(np.round(-DKY_mean,4))
        # label_h = str(t_step)
        # title = str(np.round(X,2))
        # # #file_name = str(j) + ".png"
        
        # fig = plt.figure(figsize=(10, 8))

        # ax1 = fig.add_subplot(321)
        # ax1.plot(sol[:,0], sol[:,2])
        # plt.legend([label_dky])
        # plt.xlabel('Eje x')
        # plt.ylabel('Eje z')
        # plt.title(title)
        # plt.grid(True)

        # ax2 = fig.add_subplot(322)
        # ax2.scatter(eigenvalues.real, eigenvalues.imag, marker="o")
        # plt.legend([label_eig])
        # plt.xlabel('Re')
        # plt.ylabel('Im')
        # plt.axhline(0, color="black")
        # plt.axvline(0, color="black")

        # ax1 = fig.add_subplot(323)
        # ax1.plot(t2, sol2[:,0])
        # plt.legend([label_h])
        # ax1.set_xlabel('Tiempo (s)')
        # ax1.set_ylabel('$x(t)$')

        # ax2 = fig.add_subplot(324)
        # ax2.vlines(frq2, 0, np.abs(Y2.imag))
        # plt.legend([label_potx])
        # plt.xlim(-10, 10)
        # plt.xlabel('Frecuencia (Hz)')
        # plt.ylabel('Im($Y_x$)')
        
        # ax2 = fig.add_subplot(325)
        # ax2.vlines(frq2, 0, np.abs(Y3.imag))
        # plt.legend([label_poty])
        # plt.xlim(-10, 10)
        # plt.xlabel('Frecuencia (Hz)')
        # plt.ylabel('Im($Y_y$)')
        
        # ax2 = fig.add_subplot(326)
        # ax2.vlines(frq2, 0, np.abs(Y4.imag))
        # plt.legend([label_potz])
        # plt.xlim(-10, 10)
        # plt.xlabel('Frecuencia (Hz)')
        # plt.ylabel('Im($Y_z$)')
        
        # fig.tight_layout()
        # plt.show()
        
    else: 
        DKY_mean = -1.0 
        
    # DKY es la funcion objetivo que se busca maximizar
    return DKY_mean
# ======================================================================================
# Definir el problema de optimización
class LorenzProblem(Problem):
    def __init__(self):
        super().__init__(n_var=4, n_obj=1, n_ieq_constr=0, n_eq_constr=0, xl=[0.001, 0.001, 0.001, 0.9], xu=[60, 180, 30, 0.99])
        
        # Espacios de busqueda
        # sigma = X[0]  #0.001,60
        # rho = X[1]    #0.001,180
        # beta = X[2]   #0.001,30

    def _evaluate(self, X, out, *args, **kwargs):
        out["F"] = np.array([lorenz_objective(x) for x in X])
        
        
# Guarda la salida estándar original
original_stdout = sys.stdout

class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush() # If you want the output to be visible immediately
    def flush(self) :
        for f in self.files:
            f.flush()

name = "ho4"

f = open('progress_entero_de_' + name + '.txt', 'w')
original = sys.stdout
sys.stdout = Tee(sys.stdout, f)

n_gen = 5
problem = LorenzProblem()
algorithm = DE(pop_size=10,
                initial_velocity='zero')
termination = get_termination("n_gen", n_gen)

res = minimize(problem,
               algorithm,
               termination,
               seed=1,
               save_history=True,
               verbose=True)

# Obtener las soluciones óptimas
solutions = res.X[0],res.X[1],res.X[2],res.X[3]
objectives = res.F

# Restaurar stdout
sys.stdout = original
f.close()

for i,algorithm in enumerate(res.history):
    print("Generación:", i)
    print(algorithm.pop.get("F"))
    print()   
    
    if (i == n_gen-1):
        salida = np.concatenate((algorithm.pop.get("X"),algorithm.pop.get("F")),axis=1)
        output_file = "last_gen_entero_de_" + name + ".dat"
        np.savetxt(output_file, salida, delimiter="  ")
    
# Extraer los objetivos de cada generación
objectives_history = [gen.pop.get('F').flatten() for gen in res.history]

# Convertir a un array de numpy para facilitar la escritura en un archivo
objectives_history = np.array(objectives_history)

# Guardar en un archivo de texto
np.savetxt('objectives_history_entero_de_' + name + '.txt', objectives_history)


# Imprimir las soluciones y objetivos
for solution, objective in zip(solutions, objectives):
    print("Solución:", solution)
    print("Objetivo:", objective)
    print()
    
