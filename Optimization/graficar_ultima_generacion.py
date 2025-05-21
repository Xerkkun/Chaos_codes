import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
from scipy.signal import find_peaks
from mpl_toolkits.mplot3d import axes3d

# =============================================================

def runge_kutta4(X, h, a, b, c):
    """Método de Runge Kutta de 4to orden"""
    k1 = chen(X, a, b, c)
    k2 = chen(X + h*0.5*np.array(k1), a, b, c)
    k3 = chen(X + h*0.5*np.array(k2), a, b, c)
    k4 = chen(X + h*np.array(k3), a, b, c)
    X = X + (h/6.)*(np.array(k1) + 2*np.array(k2) +
                    2*np.array(k3) + np.array(k4))

    return X

# =============================================================


def euler_forward(X, h, a, b, c):
    """Método de Euler hacia adelante"""
    dX = chen(X, a, b, c)
    X = X + h * np.array(dX)
    return X

# =============================================================


def lorenz_equations(X, sigma, rho, beta):
    """Ecuaciones del sistema de Lorenz"""
    x, y, z = X
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]
# ===============================================================================


def rossler(X, a, b, c):
    """Ecuaciones del sistema de Rossler"""
    x, y, z = X
    dx = -y - z
    dy = x + a*y
    dz = b + z*(x - c)
    return [dx, dy, dz]
# ===============================================================================

def chen(X, a, b, c):
    """Ecuaciones del sistema de Chen"""
    x, y, z = X
    dx = a*(y - x)
    dy = (c - a)*x - x*z + c*y
    dz = x*y - b*z
    return [dx, dy, dz]


# ===============================================================================


def valores_propios_lorenz(sigma, rho, beta):

    # Ajustar la h a partir de los valores propios
    x_eq, y_eq, z_eq = np.sqrt(
        np.abs(beta*(rho-1))), np.sqrt(np.abs(beta*(rho-1))), rho-1

    J = np.array(
        [[-sigma, sigma, 0], [rho - z_eq, -1, -x_eq], [y_eq, x_eq, -beta]])
    eigenvalues, _ = np.linalg.eig(J)

    vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
    vp_non_zero = vp != 0
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    # print(vp[vp_non_zero],vp_inverse)
    # t_step=0.001
    t_step = np.round(vp_min/10, 5)
    transient = int(vp_max*5/t_step)
    steady_state = int(500*vp_max)
    
    return x_eq, y_eq, z_eq, t_step, transient, eigenvalues
# -------------------------------------------------------------------------------


def valores_propios_rossler(a, b, c):
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

    vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
    vp_non_zero = vp != 0
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    # print(vp[vp_non_zero],vp_inverse)
    # t_step=0.001
    t_step = np.round(vp_min/10, 5)
    transient = int(vp_max*5/t_step)
    steady_state = int(500*vp_max)
    
    return x_eq, y_eq, z_eq, t_step, transient, eigenvalues
# -------------------------------------------------------------------------------


def valores_propios_chen(a, b, c):
    # Ajustar la h a partir de los valores propios
    x_eq = np.sqrt(2*b*c-a*b)
    y_eq = x_eq
    z_eq = np.power(x_eq, 2)/b

    J = np.array([[-a, a, 0], [c-a-z_eq, c, -x_eq], [y_eq, x_eq, -b]])
    eigenvalues, _ = np.linalg.eig(J)

    vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
    vp_non_zero = vp != 0
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    # print(vp[vp_non_zero],vp_inverse)
    # t_step=0.001
    t_step = np.round(vp_min/50, 5)
    transient = int(vp_max*5/t_step)
    steady_state = int(500*vp_max)
    
    return x_eq, y_eq, z_eq, t_step, transient, eigenvalues


# ===============================================================================

# Leer los parámetros del archivo

for semilla in range(9, 10):
    file = "last_gen_pso_chen_polaco"+str(semilla)+".dat"
    # a, b, c, DKY = np.loadtxt(file, unpack=True)
    a, b, c, DKY = np.loadtxt(file, unpack=True)

    for j in range(0, 20):
        # =============================================================
        X = np.array([a[j], b[j], c[j]])
        x_eq, y_eq, z_eq, t_step, transient, eigenvalues = valores_propios_chen(a[j], b[j], c[j])
        
        n = int(1E4)
        num_steps = transient + n

        print('Ancho de paso: ', t_step)
        print('Número de pasos: ', num_steps)
        print('Transitorio:, ', transient)

        # =============================================================
        # Condiciones iniciales (cercanas a los puntos de equilibrio)
        x0, y0, z0 = 4.246, 4.728, 13.470

        sol = np.zeros((num_steps+1, 3))
        sol[0] = [x0, y0, z0]

        # Resolver el sistema de Chen usando Euler hacia adelante
        for i in range(num_steps):
            sol[i+1] = runge_kutta4(sol[i], t_step, a[j], b[j], c[j])

            # Condición de paro por desbordamiento
            if (sol[i+1, 0] > 1E3 or sol[i+1, 0] < -1E3):
                break

        # =============================================================
        # Transformada de Fourier
        # Quitar el transitorio
        sol2 = sol[transient:, :]

        t2 = np.linspace(t_step*transient, t_step*num_steps, num_steps+1-transient)
        dt2 = t2[1] - t2[0]

        # Calcula la Transformada de Fourier
        Y2 = fft(sol2[:, 0]) / (num_steps+1-transient)
        Y3 = fft(sol2[:, 1]) / (num_steps+1-transient)
        Y4 = fft(sol2[:, 2]) / (num_steps+1-transient)
        frq2 = fftfreq(num_steps+1-transient, dt2)
        sumax = sum(abs(Y2))
        sumay = sum(abs(Y3))
        sumaz = sum(abs(Y4))

        # Encuentra los picos
        peaksx, _ = find_peaks(np.abs(Y2[0:len(sol2[:, 0])//2]))
        peaksy, _ = find_peaks(np.abs(Y3[0:len(sol2[:, 1])//2]))
        peaksz, _ = find_peaks(np.abs(Y4[0:len(sol2[:, 2])//2]))
        num_peaksx = len(peaksx)
        num_peaksy = len(peaksy)
        num_peaksz = len(peaksz)

        # Encontrar los índices de los picos que sobrepasan un valor específico
        threshold = 6
        peaks_thrx, _ = find_peaks(np.abs(Y2), height=threshold)
        peaks_thry, _ = find_peaks(np.abs(Y3), height=threshold)
        peaks_thrz, _ = find_peaks(np.abs(Y4), height=threshold)

        num_peaks_thrx = len(peaks_thrx)
        num_peaks_thry = len(peaks_thry)
        num_peaks_thrz = len(peaks_thrz)

        # Graficar
        label_potx = 'pot=' + str(np.round(sumax, 2)) + ',peaks=' + str(num_peaksx)
        label_poty = 'pot=' + str(np.round(sumay, 2)) + ',peaks=' + str(num_peaksy)
        label_potz = 'pot=' + str(np.round(sumaz, 2)) + ',peaks=' + str(num_peaksz)
        label_eig = str(np.round(eigenvalues.real, 2) +
                        np.round(eigenvalues.imag, 2) * 1j)
        label_dky = str(np.round(-DKY[j], 4))
        label_h = str(t_step)
        title = str(np.round(X, 2))
        file_name = "PSO_Chen_RK4_semilla_" + str(semilla) + "_" + str(j) + ".pdf"
        print(j)
        fig = plt.figure(figsize=(10, 8))

        ax1 = fig.add_subplot(321)
        ax1.plot(sol[:, 0], sol[:, 2], lw=0.7)
        plt.legend([label_dky])
        plt.xlabel('Eje x')
        plt.ylabel('Eje z')
        plt.title(title)

        ax2 = fig.add_subplot(322)
        ax2.scatter(eigenvalues.real, eigenvalues.imag, marker="o")
        plt.legend([label_eig])
        plt.xlabel('Re')
        plt.ylabel('Im')
        plt.axhline(0, color="black")
        plt.axvline(0, color="black")
        plt.grid(True)

        ax1 = fig.add_subplot(323)
        ax1.plot(t2, sol2[:, 0], lw=0.7)
        plt.legend([label_h])
        ax1.set_xlabel('Tiempo (s)')
        ax1.set_ylabel('$x(t)$')

        ax2 = fig.add_subplot(324)
        ax2.vlines(frq2, 0, np.abs(Y2.imag))
        plt.legend([label_potx])
        plt.xlim(-10, 10)
        plt.xlabel('Frecuencia (Hz)')
        plt.ylabel('Im($Y_x$)')

        ax2 = fig.add_subplot(325)
        ax2.vlines(frq2, 0, np.abs(Y3.imag))
        plt.legend([label_poty])
        plt.xlim(-10, 10)
        plt.xlabel('Frecuencia (Hz)')
        plt.ylabel('Im($Y_y$)')

        ax2 = fig.add_subplot(326)
        ax2.vlines(frq2, 0, np.abs(Y4.imag))
        plt.legend([label_potz])
        plt.xlim(-10, 10)
        plt.xlabel('Frecuencia (Hz)')
        plt.ylabel('Im($Y_z$)')

        fig.tight_layout()

        plt.savefig(file_name, dpi=300, bbox_inches='tight')
        plt.clf()

    # Mostrar el gráfico
    # plt.show()
