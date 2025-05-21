import numpy as np
import sys
from pymoo.core.problem import Problem
from pymoo.algorithms.soo.nonconvex.de import DE
from pymoo.termination import get_termination
from pymoo.optimize import minimize
import subprocess
from pathlib import Path
import time
from numerical_methods_v2 import runge_kutta4
from dinamics import equilibrium_points, eigenvalues, time_step, fourier_transform
from graficas import graph_optimization

start_time = time.time()
# =============================================================
# Definir la función objetivo a optimizar (Dimension de Kaplan-Yorke)
def chaos_objective(X):
    
    name = 'xu'
    
    # Parámetros del sistema 
    a, b, c, d, e, f, g = X[0], X[1], X[2], X[3], X[4], X[5], X[6]   
# =============================================================
    # Ajustar la h a partir de los valores propios
    x_eq, y_eq, z_eq = equilibrium_points(a, b, c, d, e, f, g) 
    eigenvalues_vector = eigenvalues(x_eq, y_eq, z_eq, a, b, c, d, e, f, g) 
    t_step, transient, num_steps = time_step(eigenvalues_vector)

          
    print('Ancho de paso: ', t_step)
    print('Número de pasos: ', num_steps)
    print('Transitorio:, ', transient)
# =============================================================
    # Condiciones iniciales (cercanas a los puntos de equilibrio)
    x0, y0, z0 = x_eq, -y_eq, z_eq
    
    sol = np.zeros((num_steps+1, 3))
    sol[0] = [x0, y0, z0]
    print('Condiciones iniciales: ', x0, y0, z0)
# =============================================================
    # Tiempo de integración
    t_start = 0
    t_end = num_steps*t_step
    t = np.linspace(t_start, t_end, num_steps+1)
    x_t = np.zeros((num_steps+1,1))
    x_t[:,0] = t
# =============================================================
    dd = 3 #Dimensión del sistema
    
    # Resolver el sistema de Lorenz usando Runge-Kutta 4th
    
    sol, status = runge_kutta4(sol, x_t, num_steps, t_step, a, b, c, d, e, f, g, dd, name)
    # Condición de paro por desbordamiento
    if (status == 1):
        sumax = 1
# =============================================================
    # Transformada de Fourier
    threshold = 6
    max_peak_value_x, max_peak_value_y, max_peak_value_z, num_peaksx, num_peaksy, num_peaksz, num_peaks_thrx, num_peaks_thry, num_peaks_thrz, sumax, sumay, sumaz, t2, sol2, Y2, Y3, Y4, frq2  = fourier_transform(sol, transient, t_step, num_steps, threshold)
    print("El valor más alto que alcanzan los picos para x es", max_peak_value_x)
    print("El valor más alto que alcanzan los picos para y es", max_peak_value_y)
    print("El valor más alto que alcanzan los picos para z es", max_peak_value_z)
# =============================================================
# Condición para descartar evoluciones temporales periódicas       
    # Valores propios conjugados con parte real positiva
    # Potencia de la transformada de Fourier mayor a un determinado valor
    
    #if(sumax > 150 and num_peaksx > 200):
    # TISEAN
    
    # Se guarda la evolucion temporal
    # output_file = name + ".dat"
    # np.savetxt(output_file, sol, delimiter="  ")

    k_values = (100, 100)

    DKY_sum = 0
    for i in k_values:
        # Se llama a Tisean para calcular el espectro de exponentes de lyapunov y la dimension KY
        lyap_specPATH = 'lyap_spec  ' + name + '_rk4.rnd -x' + str(transient) + ' -c2,3,4 -m3,1 -k' + str(i) +' -o' + name + '.lyaps'
        print(lyap_specPATH)
        # Ejecuta el comando de lyap_spec
        result = subprocess.run(lyap_specPATH, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Se verifica si el archivo está vacío
        with open(name + '.lyaps', "r") as file:
            file.seek(0,2)  # Mueve el cursor al final del archivo
            if file.tell() == 0:  # Verifica si el tamaño del archivo es 0 (vacío)
                DKY = 1.0
            else:
                # Si el archivo no está vacío, se obtienen los valores
                file.seek(file.tell()-9)
                DKY = float(file.read())
                file.seek(0,2)
                file.seek(file.tell()-309)
                LE = float(file.read(12))

        # Condición por si TISEAN se apendeja
        if (DKY == 3.0):
            DKY = 1.0

        DKY_sum = DKY_sum + DKY

        if(max_peak_value_z > 10.5 or num_peaks_thrz > 4):
            DKY_mean = 1.0

    DKY_mean = DKY_sum/len(k_values)

    
    file2 = open(name + "iteraciones_de.dat", "a")
    file2.write('%.5f' % X[0] + '\t' + '%.5f' % X[1] + '\t' + '%.5f' % X[2] + '\t' + '%.5f' % DKY_mean + '\n')   
    file2.close() 
    # else: DKY_mean = 1.0 
                

    print("Variables",X)
    print("Potencia", sumax)
    print("Valores propios", eigenvalues_vector)
    print("Dimension Kaplan-Yorke",DKY_mean)

    #Graficar
    graph_optimization(sumax, sumay, sumaz, num_peaksx, num_peaksy, num_peaksz, num_peaks_thrx, num_peaks_thry, num_peaks_thrz, eigenvalues_vector, DKY_mean, t_step, X, sol, t2, sol2, Y2, Y3, Y4, frq2 )

    # DKY es la funcion objetivo que se busca maximizar
    return -DKY_mean

# Definir el problema de optimización
class ChaosProblem(Problem):
    def __init__(self):
        super().__init__(n_var=7, n_obj=1, n_ieq_constr=0, n_eq_constr=0, xl=[3.001, 0.001, 0.001, 0.001, 0.001, 3.001, 0.001], xu=[4, 1, 1, 1, 1, 4, 1])
        
        # Espacios de busqueda
        # sigma = X[0]  #0.001,60
        # rho = X[1]    #0.001,180
        # beta = X[2]   #0.001,30

    def _evaluate(self, X, out, *args, **kwargs):
        out["F"] = np.array([chaos_objective(x) for x in X])
        
        
# Guarda la salida estándar original
name = 'xu'
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

f = open(name + '_progress_de.txt', 'w')
original = sys.stdout
sys.stdout = Tee(sys.stdout, f)

n_gen = 5
problem = ChaosProblem()
algorithm = DE(pop_size=10,
                initial_velocity='zero')
termination = get_termination("n_gen", n_gen)

res = minimize(problem,
               algorithm,
               termination,
               seed=5,
               save_history=True,
               verbose=True)

# Obtener las soluciones óptimas
solutions = res.X[0],res.X[1],res.X[2]
objectives = res.F

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} segundos")

# Restaurar stdout
sys.stdout = original
f.close()

for i,algorithm in enumerate(res.history):
    print("Generación:", i)
    print(algorithm.pop.get("F"))
    print()   
    
    if (i == n_gen-1):
        salida = np.concatenate((algorithm.pop.get("X"),algorithm.pop.get("F")),axis=1)
        output_file = name + "last_gen_de.dat"
        np.savetxt(output_file, salida, delimiter="  ")
    
# Extraer los objetivos de cada generación
objectives_history = [gen.pop.get('F').flatten() for gen in res.history]

# Convertir a un array de numpy para facilitar la escritura en un archivo
objectives_history = np.array(objectives_history)

# Guardar en un archivo de texto
np.savetxt(name + 'objectives_history_de.txt', objectives_history)


# Imprimir las soluciones y objetivos
for solution, objective in zip(solutions, objectives):
    print("Solución:", solution)
    print("Objetivo:", objective)
    print()
    
