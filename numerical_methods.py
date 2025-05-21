import numpy as np
from functions import lorenz_equations
#======================================================================================
def binomial_coef(alpha,mm):
    """ Cálculo de los coeficientes binomiales """
    c = np.zeros((mm+1,1))
    c[0] = 1
#    arch = open(outputfile + ".rnd","w")
    for j in range(1,mm+1):
        c[j] = c[j-1] - (c[j-1]*(1+alpha)/j)

 #       arch.write('%.15f' % c[j] + '\n')
 #   arch.close()
    return c
def grunwald_letnikov(x,h_alpha,k,alpha,x_t,nu,d,mm,decimal,sigma,rho,beta,outputfile):
    """Definición del método de Grünwald-Letnikov para la derivada de orden fraccionario"""
    # Iteraciones del método
    sum_x = np.zeros(d)
    c = binomial_coef(alpha,mm)
    arch_1 = open(outputfile + "_gl" + ".rnd","w") #"wb" para escribir archivos con formato binario
    arch_2 = open(outputfile + "_sumatoria" + ".rnd","w")
    new_x = np.zeros(d)
    aux_x = np.zeros((nu, d))
    aux_x[0,:] = x[0,:]

    for i in range(1,k+1):
        # Se calculan las sumas de los coeficientes binomiales en cada iteracion
        if nu == 1:
            for j in range(1,i+1):
                sum_x += c[j]*x[i-j,:]
            x[i,:] = lorenz_equations(x[i-1,:], sigma, rho, beta)*h_alpha - sum_x

        else:
            for j in range(1,nu+1):
                sum_x += c[j]*aux_x[nu-j,:]
            if i < mm: # Se llenan los valores de x hasta la iteracion m (maximo de memoria)
                x[i,:] = np.array(lorenz_equations(x[i-1,:], sigma, rho, beta))*h_alpha - sum_x
                aux_x[i]=x[i,:]                

           # Luego se van desechando los primeros valores de x que se habian calculado con anterioridad para ir guardando los nuevos
            # Si el  numero de iteracion (i) ya excede el valor de la memoria (m)
            # Desplazar todos los valores del vecto x una posicion hacia arriba y dejar vacio el último elemento del vector
            else:
                new_x = np.array(lorenz_equations(x[i-1,:],sigma,rho,beta))*h_alpha - sum_x
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
#======================================================================================
#Método de Runge-Kutta de 4to orden
def runge_kutta4(x, x_t, h, k, sigma, rho, beta, d, outputfile):
    """ Método de integración numérica Runge-Kutta de 4to orden """
    # Inicialización del vector de estado y tiempo
    
    arch_1 = open(outputfile + "_rk4" + ".rnd","w") #"wb" para escribir archivos con formato binario
    
    for i in range(1,k+1):
        k1 = lorenz_equations(x[i-1,:], sigma, rho, beta)
        k2 = lorenz_equations(x[i-1,:] + h*0.5*np.array(k1), sigma, rho, beta)
        k3 = lorenz_equations(x[i-1,:] + h*0.5*np.array(k2), sigma, rho, beta)
        k4 = lorenz_equations(x[i-1,:] + h*np.array(k3), sigma, rho, beta)
        x[i,:] = x[i-1,:] + (h/6.)*(np.array(k1) + 2*np.array(k2) + 2*np.array(k3) + np.array(k4))
        
        # Condición de paro por desbordamiento
        if (x[i,0]>1E3 or x[i,0]<-1E3):
            status = 0
            break
        else:
            status = 1
        if d==3:
            #if i%1000 == 0 :
            #    print(x_t[i,0],x[i,0],x[i,1],x[i,2])
            arch_1.write('%.3f' % x_t[i,0] + '\t' + '%.8f' % x[i,0] + '\t' + '%.8f' % x[i,1] + '\t' + '%.8f' % x[i,2] + '\n')
        elif d==4:
            #if i%1000 == 0 :
                #print(x_t[i,0],x[i,0],x[i,1],x[i,2],x[i,3])
            arch_1.write('%.3f' % x_t[i,0] + '\t' + '%.15f' % x[i,0] + '\t' + '%.15f' % x[i,1] + '\t' + '%.15f' % x[i,2] + '%.15f' % x[i,3] + '\n')
            
    arch_1.close()
    return x,status
#======================================================================================
def forward_euler(x,x_t, k, h, sigma, rho, beta, d, outputfile):
    
    arch_1 = open(outputfile + "_fe" + ".rnd","w") #"wb" para escribir archivos con formato binario

    for i in range(1,k+1):
        dy = lorenz_equations(x[i-1,:], sigma, rho, beta)
        x[i,:] = x[i-1,:] + h * np.array([dy])
        
        if (x[i,0]>1E3 or x[i,0]<-1E3):
            status = 0
            break
        else:
            status = 1
            
        if d==3:
            #if i%1000 == 0 :
             #   print(x_t[i,0],x[i,0],x[i,1],x[i,2])
            arch_1.write('%.3f' % x_t[i,0] + '\t' + '%.8f' % x[i,0] + '\t' + '%.8f' % x[i,1] + '\t' + '%.8f' % x[i,2] + '\n')
        elif d==4:
            #if i%1000 == 0 :
                #print(x_t[i,0],x[i,0],x[i,1],x[i,2],x[i,3])
            arch_1.write('%.3f' % x_t[i,0] + '\t' + '%.15f' % x[i,0] + '\t' + '%.15f' % x[i,1] + '\t' + '%.15f' % x[i,2] + '%.15f' % x[i,3] + '\n')
            
    arch_1.close()
    
    return x,status
#======================================================================================
#Caputo definition of fractional derivative method