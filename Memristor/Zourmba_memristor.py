import numpy as np
import math
import matplotlib.pyplot as plt 

def fx(x,a,b,fxTau):
    dx = -a*x + b*(fxTau)
    return dx

n = 10000
mu = 0.96
jj = np.linspace(0, n-1, n)
gamma_coef = np.zeros((n+1, 1))
omega_coef = np.zeros((n+1, 1))

tau = 3.5
a = 1.0
b = 4.5

gamma = 1.0
delta = 0.1

T = 100
n_tran = 2000
h = T/n

gamma_coef = np.power(jj+1,1-mu) - np.power(jj,1-mu)
omega_coef[0]=gamma_coef[0]

x = np.zeros((n+1,1))
xTau = np.zeros((n+1,1))
nTau = int(tau*n/T)

#Condición inicial
x[0] = 0.0
jx = n #Revisar que pasa cuando jx = i
for i in range(1,n+1):
    if (jx <= nTau):
        xTau[jx] = 0
        fxTau = 0.55
    else:
        xTau[jx] = x[jx - nTau]
        fxTau = gamma - delta*np.power(xTau[jx],2)

    for j in range(1,i):
        omega_coef[j] = gamma_coef[j] - gamma_coef[j-1]

    omega_coef[i] = gamma_coef[i-1]
    sum = 0
    
    for j in range(1,i):
        sum = sum + (omega_coef[j]*x[i-j])
    
    x[i] = (omega_coef[i]*x[0] - sum + np.power(h,mu)*math.gamma(2-mu)*fx(x[i-1],a,b,fxTau))/omega_coef[0]
    jx = i 

with open('resultados.txt', 'w') as file:
    file.write('j\tx\txTau\n')  # Escribir los encabezados de las columnas, separados por tabs
    for j in range(n_tran, n+1):
        # Escribir j, x[j], y xTau[j] en el archivo, separados por tabs
        file.write(f'{j}\t{x[j][0]}\t{xTau[j][0]}\n')