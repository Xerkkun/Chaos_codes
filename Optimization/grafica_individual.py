import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
from scipy.signal import find_peaks

# importar todas las funciones de pylab
from pylab import *
from mpl_toolkits.mplot3d import axes3d
from matplotlib import style
# =============================================================
#Método de Runge Kutta de 4to orden
def runge_kutta4(X, h, sigma, rho, beta):
    k1 = lorenz_equations(X, sigma, rho, beta)
    k2 = lorenz_equations(X + h*0.5*np.array(k1), sigma, rho, beta)
    k3 = lorenz_equations(X + h*0.5*np.array(k2), sigma, rho, beta)
    k4 = lorenz_equations(X + h*np.array(k3), sigma, rho, beta)
    X = X + (h/6.)*(np.array(k1) + 2*np.array(k2) + 2*np.array(k3) + np.array(k4))
    return X

# =============================================================


def euler_forward(X, h, a, b, c):
    """Método de Euler hacia adelante"""
    dX = chen(X, a, b, c)
    X = X + h * np.array(dX)
    return X

# =============================================================
# Ecuaciones del sistema de Lorenz
def lorenz_equations(X, sigma, rho, beta):
    x, y, z = X
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]
# ===============================================================================


def rossler(X, aa, bb, cc):
    x, y, z = X
    dx = -y - z
    dy = x + aa*y
    dz = bb + z*(x - cc)
    return [dx, dy, dz]
# ===============================================================================


def chen(X, aa, bb, cc):
    x, y, z = X
    dx = aa*(y - x)
    dy = (cc - aa)*x - x*z + cc*y
    dz = x*y - bb*z
    return [dx, dy, dz]


# =============================================================
aa, bb, cc = 28.21888891274968714,  0.7066636808327453334,  24.75217691457946856
DKY = 2.3138
# =============================================================
X = np.array([aa,bb,cc])
# print(X)
# Ajustar la h a partir de los valores propios
# Ajustar la h a partir de los valores propios
x_eq = np.sqrt(2*bb*cc-aa*bb)
y_eq = x_eq
z_eq = np.power(x_eq, 2)/bb

J = np.array([[-aa, aa, 0], [cc-aa-z_eq, cc, -x_eq], [y_eq, x_eq, -bb]])
eigenvalues, _ = np.linalg.eig(J)

vp = np.array([eigenvalues.real, eigenvalues.imag], dtype=float)
vp_non_zero = vp != 0
vp_inverse = 1/np.abs(vp[vp_non_zero])
vp_min = np.amin(vp_inverse)
vp_max = np.amax(vp_inverse)
# print(vp[vp_non_zero],vp_inverse)
# t_step=0.001
t_step = np.round(vp_min/50, 5)
n = int(1E4)

transient=int(vp_max*5/t_step)
num_steps=transient+n

# print('Ancho de paso: ', t_step)
# print('Número de pasos: ', num_steps)
# print('Transitorio:, ', transient)
# =============================================================
#Condiciones iniciales (cercanas a los puntos de equilibrio)
x0, y0, z0 = 4.246, 4.728, 13.470

sol = np.zeros((num_steps+1, 3))
sol[0] = [x0, y0, z0]

# Resolver el sistema de Lorenz usando Runge-Kutta 4th
for i in range(num_steps):
    sol[i+1] = euler_forward(sol[i],t_step,aa,bb,cc)
    
    # Condición de paro por desbordamiento
    if (sol[i+1,0]>1E3 or sol[i+1,0]<-1E3):
        suma = 1
        break
# =============================================================
# Transformada de Fourier
# Quitar el transitorio
sol2 = sol[transient:,:]

t2 = np.linspace(t_step*transient, t_step*num_steps, num_steps+1-transient)
dt2 = t2[1] - t2[0]

# Calcula la Transformada de Fourier
Y2 = fft(sol2[:,0]) / (num_steps+1-transient)  # Transformada normalizada
Y3 = fft(sol2[:,1]) / (num_steps+1-transient)  # Transformada normalizada
Y4 = fft(sol2[:,2]) / (num_steps+1-transient)  # Transformada normalizada
frq2 = fftfreq(num_steps+1-transient, dt2) 
sumax = sum(abs(Y2))    
sumay = sum(abs(Y3))    
sumaz = sum(abs(Y4))    

# Encuentra los picos
peaksx, _ = find_peaks(np.abs(Y2[0:len(sol2[:,0])//2]))
peaksy, _ = find_peaks(np.abs(Y3[0:len(sol2[:,1])//2]))
peaksz, _ = find_peaks(np.abs(Y4[0:len(sol2[:,2])//2]))
num_peaksx = len(peaksx)
num_peaksy = len(peaksy)
num_peaksz = len(peaksz)

# Encuentra picos a partir de un umbral
# Calcular la magnitud de la transformada de Fourier
magnitudex = np.abs(Y2)
magnitudey = np.abs(Y3)
magnitudez = np.abs(Y4)

# Encontrar los índices de los picos que sobrepasan un valor específico
threshold = 6
peaks_thrx, _ = find_peaks(magnitudex, height=threshold)
peaks_thry, _ = find_peaks(magnitudey, height=threshold)
peaks_thrz, _ = find_peaks(magnitudez, height=threshold)
#    print("Indices de los picos que sobrepasan el umbral:", peaks)    

num_peaks_thrx = len(peaks_thrx)

num_peaks_thry = len(peaks_thry)
num_peaks_thrz = len(peaks_thrz)    

# Calcular el valor más alto que alcanzan los picos
# max_peak_value_x = np.max(magnitudex[peaksx])
# max_peak_value_y = np.max(magnitudey[peaksy])
# max_peak_value_z = np.max(magnitudez[peaksz])

# print("El valor más alto que alcanzan los picos para x es", max_peak_value_x)
# print("El valor más alto que alcanzan los picos para y es", max_peak_value_y)
# print("El valor más alto que alcanzan los picos para z es", max_peak_value_z)

# =============================================================
# Graficar 

# subplot(311)
# p1,=plot(t2,sol2[:,0],"m",lw=0.7)
# ylabel("x")
# plt.xticks([])

# subplot(312)
# p2,=plot(t2,sol2[:,1],"m",lw=0.7)
# ylabel("y")
# plt.xticks([])

# subplot(313)
# p3,=plot(t2,sol2[:,2],"m",lw=0.7)
# ylabel("z")
# xlabel("Tiempo")
# plt.savefig("time1.pdf",dpi=300,bbox_inches= 'tight')

# clf()

# p1,=plot(sol2[:,0],sol2[:,1],"m",lw=0.3)
# xlabel("x")
# ylabel("y")
# plt.savefig("x-y1.pdf",dpi=300,bbox_inches= 'tight')
# clf()

# p2,=plot(sol2[:,0],sol2[:,2],"m",lw=0.3)
# xlabel("x")
# ylabel("z")
# plt.savefig("x-z1.pdf",dpi=300,bbox_inches= 'tight')
# clf()

# p3,=plot(sol2[:,1],sol2[:,2],"m",lw=0.3)
# xlabel("y")
# ylabel("z")
# plt.savefig("y-z1.pdf",dpi=300,bbox_inches= 'tight')
# clf()

# Atractor 3d
# fig = plt.figure()
# ax1 = fig.add_subplot(111,projection='3d')
# ax1.plot(sol2[:,0],sol2[:,1],sol2[:,2],"m",lw=0.2)
# ax1.set_xlabel('x')
# ax1.set_ylabel('y')
# ax1.set_zlabel('z')
# plt.savefig("atractor3d1.pdf",dpi=300,bbox_inches= 'tight')
# clf()
    
# label_eig = str(np.round(eigenvalues.real, 2) + np.round(eigenvalues.imag, 2) * 1j)
# file_name = "example2.pdf"

# fig = plt.figure(figsize=(10, 8))

# ax1 = fig.add_subplot(221)
# ax1.plot(sol[:,0], sol[:,2],lw=0.7)
# plt.xlabel('x')
# plt.ylabel('z')

# ax2 = fig.add_subplot(222)
# ax2.scatter(eigenvalues.real, eigenvalues.imag, marker="o")
# plt.xlabel('Re')
# plt.ylabel('Im')
# plt.axhline(0, color="black")
# plt.axvline(0, color="black")
# plt.xlim(-50,2)
# plt.grid(True)

# plt.savefig(file_name,dpi=300,bbox_inches='tight')
# plt.clf()

# file_name = "example3.pdf"

# fig = plt.figure(figsize=(10, 5))

# ax1 = fig.add_subplot(221)
# ax1.plot(t2, sol2[:,0],lw=0.7)
# ax1.set_xlabel('Iterations')
# ax1.set_ylabel('$x(t)$')

# ax2 = fig.add_subplot(222)
# ax2.vlines(frq2, 0, np.abs(Y2.imag))
# plt.xlim(-10, 10)
# # plt.ylim(0,max_peak_value_x)
# plt.xlabel('Frecuency')
# plt.ylabel('X')

# fig.tight_layout()

# plt.savefig(file_name,dpi=300,bbox_inches='tight')
# plt.clf()

# =============================================================
    # Graficar 
label_eig = str(np.round(eigenvalues.real, 2) + np.round(eigenvalues.imag, 2) * 1j)
label_dky = str(np.round(DKY,4))
label_h = 'h = '+str(t_step)
file_name = "best_chen.pdf"

fig = plt.figure(figsize=(10, 8))

ax1 = fig.add_subplot(221)
ax1.plot(sol[:,0], sol[:,2],lw=0.7)
plt.legend([label_dky])
plt.xlabel('x')
plt.ylabel('z')

ax2 = fig.add_subplot(222)
ax2.scatter(eigenvalues.real, eigenvalues.imag, marker="o")
plt.legend([label_eig])
plt.xlabel('Re')
plt.ylabel('Im')
plt.axhline(0, color="black")
plt.axvline(0, color="black")
plt.grid(True)

ax1 = fig.add_subplot(223)
ax1.plot(t2, sol2[:,0],lw=0.7)
plt.legend([label_h])
ax1.set_xlabel('Iterations')
ax1.set_ylabel('$x(t)$')

ax2 = fig.add_subplot(224)
ax2.vlines(frq2, 0, np.abs(Y2.imag))
plt.xlim(-5, 5)
# plt.ylim(0,max_peak_value_x)
plt.xlabel('Frecuency (Hz)')
plt.ylabel('X')

fig.tight_layout()

plt.savefig(file_name,dpi=300,bbox_inches='tight')
plt.clf()