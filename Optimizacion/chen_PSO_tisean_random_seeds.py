"""Date: 06/03/2024 
    TISEAN Lyapunov's espectrum &
    Runs with random seed"""

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from pymoo.core.problem import Problem
from pymoo.algorithms.soo.nonconvex.pso import PSO
from pymoo.termination import get_termination
from pymoo.optimize import minimize
import subprocess
from pathlib import Path
from scipy.fftpack import fft, fftfreq
from scipy.signal import find_peaks
import time

# global_i = 0
start_time = time.time()
# =============================================================
# Método de Runge Kutta de 4to orden


def runge_kutta4(X, h, a, b, c):
    k1 = chen(X, a, b, c)
    k2 = chen(X + h*0.5*np.array(k1), a, b, c)
    k3 = chen(X + h*0.5*np.array(k2), a, b, c)
    k4 = chen(X + h*np.array(k3), a, b, c)
    X = X + (h/6.)*(np.array(k1) + 2*np.array(k2) +
                    2*np.array(k3) + np.array(k4))
    return X
# =============================================================
# Ecuaciones del sistema de Lorenz


def chen(X, a, b, c):
    x, y, z = X
    dx = a*(y - x)
    dy = (c - a)*x - x*z + c*y
    dz = x*y - b*z
    return [dx, dy, dz]
# =============================================================
# Definir la función objetivo a optimizar (Dimension de Kaplan-Yorke)


def lorenz_objective(X):
    # global global_i

    # Si son muchos valores, decimar para ingresar a tisean
    # Parámetros del sistema de Lorenz
    a = X[0]  # 0,60
    b = X[1]  # 0.001,180
    c = X[2]  # 0.001,30
# =============================================================
   # Ajustar la h a partir de los valores propios
    if (b*(2*c-a) >= 0):
        x_eq = np.sqrt(2*b*c-a*b)
        y_eq = x_eq
        z_eq = np.power(x_eq, 2)/b
    else:
        x_eq = 0
        y_eq = 0
        z_eq = 0

    J = np.array([[-a, a, 0], [c-a-z_eq, c, -x_eq], [y_eq, x_eq, -b]])
    eigenvalues, _ = np.linalg.eig(J)

    vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
    vp_non_zero = vp != 0
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    # print(vp[vp_non_zero],vp_inverse)
    # t_step=0.001
    t_step = np.round(vp_min/10, 5)

    transient = int(vp_max*10/t_step)
    steady_state = int(5000*vp_max)  # int(1E5)
    num_steps = transient+steady_state
    n = int(1E4)
    # print('Ancho de paso: ', t_step)
    # print('Número de pasos: ', num_steps)
    # print('Transitorio:, ', transient)
# =============================================================
    # Condiciones iniciales (cercanas a los puntos de equilibrio)
    x0, y0, z0 = 4.246, 4.728, 13.470
# =============================================================
    # Tiempo de integración
    # t_start = 0
    # t_end = 50
    # t_step = 0.01
    # num_steps = int((t_end - t_start) / t_step)

    sol = np.zeros((num_steps+1, 3))
    sol[0] = [x0, y0, z0]

    # Resolver el sistema de Lorenz usando Forward Euler
    # for i in range(num_steps):
    #     dx, dy, dz = chen(sol[i], a, b, c)
    #     sol[i+1] = sol[i] + t_step * np.array([dx, dy, dz])

    # Resolver el sistema de Lorenz usando Runge-Kutta 4th
    for i in range(num_steps):
        sol[i+1] = runge_kutta4(sol[i], t_step, a, b, c)

        # Condición de paro por desbordamiento
        if (sol[i+1, 0] > 1E3 or sol[i+1, 0] < -1E3):
            sumax = 1
            break
# =============================================================
    # Transformada de Fourier
    # Quitar el transitorio
    sol2 = sol[transient:, :]
    # t2 = np.linspace(t_step*transient, t_step*num_steps, num_steps+1-transient)

    # dt2 = t2[1] - t2[0]

    Y2 = fft(sol2[:, 0]) / (num_steps+1-transient)  # Transformada normalizada
    # Y3 = fft(sol2[:,1]) / (num_steps+1-transient)  # Transformada normalizada
    Y4 = fft(sol2[:, 2]) / (num_steps+1-transient)  # Transformada normalizada
    # frq2 = fftfreq(num_steps+1-transient, dt2)
    sumax = sum(abs(Y2))
    # sumay = sum(abs(Y3))
    sumaz = sum(abs(Y4))

    # Encuentra los picos
    peaksx, _ = find_peaks(np.abs(Y2[0:len(sol2[:, 0])//2]))
    # peaksy, _ = find_peaks(np.abs(Y3[0:len(sol2[:,1])//2]))
    peaksz, _ = find_peaks(np.abs(Y4[0:len(sol2[:, 2])//2]))
    num_peaksx = len(peaksx)
    # num_peaksy = len(peaksy)
    # num_peaksz = len(peaksz)

    # Encuentra picos a partir de un umbral
    # Calcular la magnitud de la transformada de Fourier
    # magnitudex = np.abs(Y2)
    # magnitudey = np.abs(Y3)
    magnitudez = np.abs(Y4)

    # Encontrar los índices de los picos que sobrepasan un valor específico
    threshold = 6
    # peaks_thrx, _ = find_peaks(magnitudex, height=threshold)
    # peaks_thry, _ = find_peaks(magnitudey, height=threshold)
    peaks_thrz, _ = find_peaks(magnitudez, height=threshold)

    # num_peaks_thrx = len(peaks_thrx)
    # num_peaks_thry = len(peaks_thry)
    num_peaks_thrz = len(peaks_thrz)

   # Calcular el valor más alto que alcanzan los picos
    # max_peak_value_x = np.max(magnitudex[peaksx])
    # max_peak_value_y = np.max(magnitudey[peaksy])
    #

    # print("El valor más alto que alcanzan los picos para x es", max_peak_value_x)
    # print("El valor más alto que alcanzan los picos para y es", max_peak_value_y)
    # print("El valor más alto que alcanzan los picos para z es", max_peak_value_z)

# =============================================================
# Condición para descartar evoluciones temporales periódicas
    # Valores propios conjugados con parte real positiva
    # Potencia de la transformada de Fourier mayor a un determinado valor

    if (sumax > 20 and sumaz < 100):
        # TISEAN

        # Se guarda la evolucion temporal
        output_file = "chen_pso.dat"
        np.savetxt(output_file, sol, delimiter="  ")

        # print("Datos guardados en", output_file)

        lyap_specPATH = 'lyap_spec chen_pso.dat -x' + \
            str(transient) + ' -c1,2,3 -m3,1 -k100 -r0.01 -osalida_pso_chen.lyaps'

        # check_output(str(lyap_specPATH), shell=True)

        # Ejecuta el comando de lyap_spec y redirige la salida estándar y la salida de error a /dev/null
        result = subprocess.run(
            lyap_specPATH, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
       
        archivo = "salida_pso_chen.lyaps"
        
        # Verificar si el archivo está vacío
        if os.path.getsize(archivo) > 0:
            with open(archivo, "r") as file:
                for linea in file:
                    partes = linea.strip().split()  # Dividir la línea en partes
                    if len(partes) >= 3:  # Asegurar que la línea tiene al menos 3 partes
                        # Obtener el valor de la tercera columna
                        valor = partes[2]
                        if valor.lower() == "-nan":  # Verificar si el valor es -nan
                            DKY = 1.0
                        else:
                            # Mover el cursor al final del archivo
                            file.seek(0, 2)
                            if file.tell() >= 9:  # Verificar si el archivo tiene al menos 9 bytes para leer
                                file.seek(file.tell()-9)
                                DKY = float(file.read())

        else:
            print("El archivo está vacío.")
            # Aquí puedes manejar el caso de archivo vacío, por ejemplo, asignando un valor predeterminado a DKY
            DKY = 1.0  # O asignar un valor predeterminado adecuado

        if DKY >= 2.4:
            DKY = 1.0
        elif DKY == 3.0:
            DKY = 1.0

    else:
        DKY = 1.0

    # print("Variables",X)
    # print("Potencia", suma)
    # print("Valores propios", eigenvalues)
    # print("Dimension Kaplan-Yorke",-DKY_mean)

    # Graficar
    # label_potx = 'pot=' + str(np.round(sumax,2)) + ',peaks=' + str(num_peaksx) + ',' + str(num_peaks_thrx)
    # label_poty = 'pot=' + str(np.round(sumay,2)) + ',peaks=' + str(num_peaksy) + ',' + str(num_peaks_thry)
    # label_potz = 'pot=' + str(np.round(sumaz,2)) + ',peaks=' + str(num_peaksz) + ',' + str(num_peaks_thrz)
    # label_eig = str(np.round(eigenvalues.real, 2) + np.round(eigenvalues.imag, 2) * 1j)
    # label_dky = str(np.round(-DKY_mean,4))
    # label_h = str(t_step)
    # title = str(np.round(X,2))
    # #file_name = str(j) + ".png"

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

    # DKY es la funcion objetivo que se busca maximizar
    return -DKY

# Definir el problema de optimización


class LorenzProblem(Problem):
    def __init__(self):
        super().__init__(n_var=3, n_obj=1, n_ieq_constr=0,
                         n_eq_constr=0, xl=[20.0001, 0.0001, 20.0001], xu=[45, 20, 40])

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
            f.flush()  # If you want the output to be visible immediately

    def flush(self):
        for f in self.files:
            f.flush()


# Aquí va tu código
n_gen = 20
problem = LorenzProblem()
algorithm = PSO(pop_size=20,
                initial_velocity='zero',
                w=0.4,
                c1=1.8,
                c2=2.2)
termination = get_termination("n_gen", n_gen)

random_seed = np.random.randint(1, 20, size=10)

for run in range(8, 10):
    f = open("progress_pso_chen_tisean" + str(run) + ".txt", 'w')
    original = sys.stdout
    sys.stdout = Tee(sys.stdout, f)

    res = minimize(problem,
                   algorithm,
                   termination,
                   seed=random_seed[run],
                   save_history=True,
                   verbose=True)

    # Obtener las soluciones óptimas
    solutions = res.X[0], res.X[1], res.X[2]
    objectives = res.F

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} segundos")

    # Restaurar stdout
    sys.stdout = original
    f.close()

    last_generation = res.history[-1]
    solutions_last_gen = last_generation.pop.get('X')

    salida = np.concatenate(
        (solutions_last_gen, last_generation.pop.get("F")), axis=1)
    output_file = "last_gen_pso_chen_tisean" + str(run) + ".dat"
    np.savetxt(output_file, salida, delimiter="  ")
