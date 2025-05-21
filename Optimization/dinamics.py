import numpy as np
from scipy.fftpack import fft, fftfreq
from scipy.signal import find_peaks

def equilibrium_points(a, b, c, d, e, f, g):
# =============================================================
    # Ajustar la h a partir de los valores propios
    x_eq = (-a*e + np.sqrt((a*e)**2 + 4*(b*e + c*d)*(c*f)))/(2*(b*e + c*d))
    y_eq, z_eq = (a*x_eq + b*x_eq**2), 0
    return x_eq, y_eq, z_eq

def eigenvalues(x_eq, y_eq, z_eq, a, b, c, d, e, f, g):
    J = np.array([ [0, 0, g], [2*d*x_eq, 2*e*y_eq, 0], [-a-2*b*x_eq, 2*c*y_eq, 0]  ])
    eigenvalues, _ = np.linalg.eig(J)
    
    return eigenvalues
    
def time_step(eigenvalues):
    vp = np.array( [eigenvalues.real,eigenvalues.imag], dtype=float )
    vp_non_zero = vp != 0 
    vp_inverse = 1/np.abs(vp[vp_non_zero])
    vp_min = np.amin(vp_inverse)
    vp_max = np.amax(vp_inverse)
    #print(vp[vp_non_zero],vp_inverse)

    t_step=np.round(vp_min/10,5)
    
    #transient=int(vp_max*5/t_step)
    #transient=int(vp_max)
    transient = int(5E3)
    steady_state = int(1E4)
    num_steps=transient+steady_state
    return t_step, transient, num_steps

def fourier_transform(sol, transient, t_step, num_steps, threshold):
    sol2 = sol[transient:,:]
    t2 = np.linspace(t_step*transient, t_step*num_steps, num_steps+1-transient)

    dt2 = t2[1] - t2[0]

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
    peaks_thrx, _ = find_peaks(magnitudex, height=threshold)
    peaks_thry, _ = find_peaks(magnitudey, height=threshold)
    peaks_thrz, _ = find_peaks(magnitudez, height=threshold)
    
    num_peaks_thrx = len(peaks_thrx)
    num_peaks_thry = len(peaks_thry)
    num_peaks_thrz = len(peaks_thrz)
    
   #Calcular el valor más alto que alcanzan los picos
    max_peak_value_x = np.max(magnitudex[peaksx])
    max_peak_value_y = np.max(magnitudey[peaksy])
    max_peak_value_z = np.max(magnitudez[peaksz])

    return max_peak_value_x, max_peak_value_y, max_peak_value_z, num_peaksx, num_peaksy, num_peaksz, num_peaks_thrx, num_peaks_thry, num_peaks_thrz, sumax, sumay, sumaz, t2, sol2, Y2, Y3, Y4, frq2 