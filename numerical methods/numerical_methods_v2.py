import numpy as np
from functions import rossler

def write_to_file(filename, data):
    """Write data to a file."""
    with open(filename, "w") as file:
        for line in data:
            file.write(line + '\n')

def binomial_coef(alpha, mm):
    """Calculate binomial coefficients."""
    c = np.zeros((mm+1,1))
    c[0] = 1
    for j in range(1, mm+1):
        c[j] = c[j-1] - (c[j-1] * (1 + alpha) / j)
    return c

def grunwald_letnikov(x, h_alpha, k, alpha, x_t, nu, dd, mm, a, b, c, d, outputfile):
    """Grunwald-Letnikov method for fractional order derivative."""    
    sum_x = np.zeros(dd)
    c = binomial_coef(alpha,mm)
    new_x = np.zeros(dd)
    aux_x = np.zeros((nu, dd))
    aux_x[0,:] = x[0,:]
    data_1, data_2 = [], []

    for i in range(1,k+1):
        # Se calculan las sumas de los coeficientes binomiales en cada iteracion
        if nu == 1:
            for j in range(1,i+1):
                sum_x += c[j]*x[i-j,:]
            x[i,:] = rossler(x[i-1,:], a, b, c)*h_alpha - sum_x

        else:
            for j in range(1,nu+1):
                sum_x += c[j]*aux_x[nu-j,:]
            if i < mm: # Se llenan los valores de x hasta la iteracion m (maximo de memoria)
                x[i,:] = np.array(rossler(x[i-1,:], a, b, c))*h_alpha - sum_x
                aux_x[i]=x[i,:]                

           # Luego se van desechando los primeros valores de x que se habian calculado con anterioridad para ir guardando los nuevos
            # Si el  numero de iteracion (i) ya excede el valor de la memoria (m)
            # Desplazar todos los valores del vecto x una posicion hacia arriba y dejar vacio el último elemento del vector
            else:
                new_x = np.array(rossler(x[i-1,:], a, b, c))*h_alpha - sum_x
                aux_x = np.concatenate((aux_x[1:],[new_x]),axis=0)
                x[i,:] = new_x        
                
        # Condición de paro por desbordamiento
        if (x[i,0]>1E3 or x[i,0]<-1E3):
            status = 0
            print("saludos")
            break
        else:
            status = 1
 
        if dd == 3:
            data_1.append(f'{x_t[i,0]:.3f}\t{x[i,0]:.8f}\t{x[i,1]:.8f}\t{x[i,2]:.8f}')
            data_2.append(f'{sum_x[0]:.8f}\t{sum_x[1]:.8f}\t{sum_x[2]:.8f}')
        elif dd == 4:
            data_1.append(f'{x_t[i,0]:.3f}\t{x[i,0]:.15f}\t{x[i,1]:.15f}\t{x[i,2]:.15f}\t{x[i,3]:.15f}')
            data_2.append(f'{sum_x[0]:.15f}\t{sum_x[1]:.15f}\t{sum_x[2]:.15f}\t{sum_x[3]:.15f}')

        sum_x = np.zeros(d)
        
        write_to_file(outputfile + "_gl.rnd", data_1)
        write_to_file(outputfile + "_sumatoria.rnd", data_2)
        
    return x, status

def runge_kutta4(x, x_t, k, h, a, b, c, dd, outputfile):
    """4th order Runge-Kutta numerical integration method."""
    data = []
    
    for i in range(1, k+1):
        k1 = rossler(x[i-1,:], a, b, c)
        k2 = rossler(x[i-1,:] + h*0.5*np.array(k1), a, b, c)
        k3 = rossler(x[i-1,:] + h*0.5*np.array(k2), a, b, c)
        k4 = rossler(x[i-1,:] + h*np.array(k3), a, b, c)
        x[i,:] = x[i-1,:] + (h/6.)*(np.array(k1) + 2*np.array(k2) + 2*np.array(k3) + np.array(k4))
        
        # Condición de paro por desbordamiento
        if (x[i,0]>1E3 or x[i,0]<-1E3):
            status = 0
            break
        else:
            status = 1
        
        if dd == 3:
            data.append(f'{x_t[i,0]:.3f}\t{x[i,0]:.8f}\t{x[i,1]:.8f}\t{x[i,2]:.8f}')
        elif dd == 4:
            data.append(f'{x_t[i,0]:.3f}\t{x[i,0]:.15f}\t{x[i,1]:.15f}\t{x[i,2]:.15f}\t{x[i,3]:.15f}')
    
    write_to_file(outputfile + "_rk4.rnd", data)
    
    return x, status

def forward_euler(x, x_t, k, h, a, b, c, dd, outputfile):
    """Forward Euler method."""
    data = []
    
    for i in range(1, k+1):
        dy = rossler(x[i-1,:], a, b, c)
        x[i,:] = x[i-1,:] + h * np.array([dy])
        
        if (x[i,0]>1E3 or x[i,0]<-1E3):
            status = 0
            break
        else:
            status = 1
            
        
        if dd == 3:
            data.append(f'{x_t[i,0]:.3f}\t{x[i,0]:.8f}\t{x[i,1]:.8f}\t{x[i,2]:.8f}')
        elif dd == 4:
            data.append(f'{x_t[i,0]:.3f}\t{x[i,0]:.15f}\t{x[i,1]:.15f}\t{x[i,2]:.15f}\t{x[i,3]:.15f}')
    
    write_to_file(outputfile + "_fe.rnd", data)
    
    return x, status

#EN PROCESO
def caputo_derivative(x, x_t, k, h, a, b, c, dd, outputfile):
    """Caputo definition of fractional derivative method."""
    data = []
    
    for i in range(1, k+1):
        
        if (x[i,0]>1E3 or x[i,0]<-1E3):
            status = 0
            break
        else:
            status = 1
            
        
        if dd == 3:
            data.append(f'{x_t[i,0]:.3f}\t{x[i,0]:.8f}\t{x[i,1]:.8f}\t{x[i,2]:.8f}')
        elif dd == 4:
            data.append(f'{x_t[i,0]:.3f}\t{x[i,0]:.15f}\t{x[i,1]:.15f}\t{x[i,2]:.15f}\t{x[i,3]:.15f}')
    
    write_to_file(outputfile + "_ca.rnd", data)
        
    return x, status