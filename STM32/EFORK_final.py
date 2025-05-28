import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Para gráficas 3D

# Seleccionar tipo de dato de mayor precisión disponible
dtype = np.float128 if hasattr(np, 'float128') else np.float64

# Definición de los sistemas caóticos
def lorenz_system(x, y, z):
    sigma, rho, beta = 10.0, 28.0, 8.0/3.0
    return np.array([sigma*(y-x), x*(rho-z)-y, x*y-beta*z])

def rossler_system(x, y, z):
    a, b, c = 0.2, 0.2, 5.7
    return np.array([-y-z, x+a*y, b+z*(x-c)])

def chen_system(x, y, z):
    u, v, w = 7.5, 1.0, 5.0
    return np.array([u*(y-x), (w-u)*x-x*z+w*y, x*y-v*z])

# Término de memoria fraccional vectorizado para EFORK
def memory_fractional(k, t, arr, vtn, h, alpha, nu):
    start = max(0, k-nu)
    idx = np.arange(start, k)
    t0, t1 = vtn[idx], vtn[idx+1]
    v1 = (t - t0)**(1-alpha)
    v2 = (t - t1)**(1-alpha)
    coeff = (arr[idx+1] - arr[idx])/(h*gamma(2-alpha))
    return np.sum(coeff * (v1 - v2))

# Función genérica EFORK con término de memoria
def fn_efork(k, t, x, y, z, vtn, vxn, vyn, vzn, h, alpha, nu, system_func):
    dxdt, dydt, dzdt = system_func(x, y, z)
    mem_x = memory_fractional(k, t, vxn, vtn, h, alpha, nu)
    mem_y = memory_fractional(k, t, vyn, vtn, h, alpha, nu)
    mem_z = memory_fractional(k, t, vzn, vtn, h, alpha, nu)
    return np.array([dxdt-mem_x, dydt-mem_y, dzdt-mem_z])

# Menú y simulación principal usando inputs con dtype
def main():
    # Selección de sistema EFORK
    dic = {'1':'Lorenz','2':'Rossler','3':'Chen'}
    print("Seleccione sistema EFORK a simular:")
    for k,v in dic.items(): print(f"  {k}: {v}")
    choice = input("Opción (1/2/3): ")
    system_name = dic.get(choice)
    if not system_name:
        print("Selección inválida.")
        return
    system_func = {'Lorenz':lorenz_system,'Rossler':rossler_system,'Chen':chen_system}[system_name]

    # Solicitar parámetros al usuario
    alpha = dtype(float(input("Ingrese orden fraccionario alpha: ")))
    h     = dtype(float(input("Ingrese paso h: ")))
    Lm    = dtype(float(input("Ingrese longitud de memoria Lm: ")))
    t_f   = dtype(float(input("Ingrese tiempo final t_f: ")))
    x0_str = input("Ingrese condiciones iniciales x0,y0,z0 separadas por coma: ")
    x0_arr = np.array([float(val) for val in x0_str.split(',')], dtype=dtype)
    x0, y0, z0 = x0_arr

    # Parámetros de EFORK derivados
    dt = h
    nu = int(Lm / h)
    ha = h**alpha           # h^alpha
    # Constantes EFORK necesarias
    c2 = (1.0/(2*gamma(1+alpha)))**(1.0/alpha)
    c3 = (1.0/(4*gamma(1+alpha)))**(1.0/alpha)

    a21 = 1.0/(2*gamma(alpha+1)**2)
    a31 = ( gamma(alpha+1)**2 * gamma(2*alpha+1) + 2*gamma(2*alpha+1)**2 - gamma(3*alpha+1) ) / ( 4*gamma(alpha+1)**2 * (2*gamma(2*alpha+1)**2 - gamma(3*alpha+1)) )
    a32 = - gamma(2*alpha+1) / ( 4*(2*gamma(2*alpha+1)**2 - gamma(3*alpha+1)) )

    w1 = (8.0*gamma(1+alpha)**3 * gamma(1+2*alpha)**2 - 6.0*gamma(1+alpha)**3 * gamma(1+3*alpha) + gamma(1+2*alpha)*gamma(1+3*alpha)) / (gamma(1+alpha)*gamma(1+2*alpha)*gamma(1+3*alpha))
    w2 = 2.0*gamma(1+alpha)**2 * (4.0*gamma(1+2*alpha)**2 - gamma(1+3*alpha)) / (gamma(1+2*alpha)*gamma(1+3*alpha))
    w3 = -8.0*gamma(1+alpha)**2 * (2.0*gamma(1+2*alpha)**2 - gamma(1+3*alpha)) / (gamma(1+2*alpha)*gamma(1+3*alpha))

    # Convertir tiempo final a número de pasos
    dt = h
    N1 = int(np.ceil(t_f / h))
    # Inicializar arrays de memoria
    vtn = np.zeros(N1+1)
    vxn = np.zeros(N1+1)
    vyn = np.zeros(N1+1)
    vzn = np.zeros(N1+1)
    vtn[0], vxn[0], vyn[0], vzn[0] = 0.0, x0, y0, z0
    n, x_n, y_n, z_n = 0, x0, y0, z0
    # Primer paso sin término de memoria
    K1 = ha * system_func(x_n, y_n, z_n)
    K2 = ha * system_func(x_n + a21*K1[0], 
                          y_n + a21*K1[1], 
                          z_n + a21*K1[2])
    K3 = ha * system_func(x_n + a31*K2[0] + a32*K1[0],
                          y_n + a31*K2[1] + a32*K1[1],
                          z_n + a31*K2[2] + a32*K1[2])
    
    x_n1 = x_n + w1*K1[0] + w2*K2[0] + w3*K3[0]
    y_n1 = y_n + w1*K1[1] + w2*K2[1] + w3*K2[1] 
    z_n1 = z_n + w1*K1[2] + w2*K2[2] + w3*K3[2]
    
    
    tn1 = (n+1)*h
    vtn[1], vxn[1], vyn[1], vzn[1] = tn1, x_n1, y_n1, z_n1
    x_n, y_n, z_n = x_n1, y_n1, z_n1

    # Bucle principal de EFORK
    for n in range(1, N1):
        tn = n * h

        K1 = ha * fn_efork(n, tn, x_n, y_n, z_n, vtn, vxn, vyn, vzn, h, alpha, nu, system_func)
        K2 = ha * fn_efork(n, tn + c2*h,
                           x_n + a21*K1[0], 
                           y_n + a21*K1[1], 
                           z_n + a21*K1[2],
                           vtn, vxn, vyn, vzn, h, alpha, nu, system_func)
        K3 = ha * fn_efork(n, tn + c3*h,
                           x_n + a31*K2[0] + a32*K1[0],
                           y_n + a31*K2[1] + a32*K1[1],
                           z_n + a31*K2[2] + a32*K1[2],
                           vtn, vxn, vyn, vzn, h, alpha, nu, system_func)
        
        x_n1 = x_n + w1*K1[0] + w2*K2[0] + w3*K3[0]
        y_n1 = y_n + w1*K1[1] + w2*K2[1] + w3*K3[1]
        z_n1 = z_n + w1*K1[2] + w2*K2[2] + w3*K3[2]

        tn1 = (n+1)*h
        vtn[n+1], vxn[n+1], vyn[n+1], vzn[n+1] = tn1, x_n1, y_n1, z_n1

        x_n, y_n, z_n = x_n1, y_n1, z_n1

    # Guardar resultados
    datos = np.column_stack((vtn, vxn, vyn, vzn))
    out_txt = f"EFORK_{system_name}_alpha{alpha}.txt"
    np.savetxt(out_txt, datos, delimiter='\t', header='t\tx\ty\tz', comments='')
    print(f"Resultados guardados en: {out_txt}")

    pdf2d = f"{system_name}_atractores_EFORK_alpha{alpha}.pdf"
    plt.figure(figsize=(12, 4))
    # Proyección xy
    plt.subplot(1, 3, 1)
    plt.plot(vxn, vyn, 'm', lw=0.3)
    plt.xlabel('x'); plt.ylabel('y')
    # Proyección xz
    plt.subplot(1, 3, 2)
    plt.plot(vxn, vzn, 'm', lw=0.3)
    plt.xlabel('x'); plt.ylabel('z')
    # Proyección yz
    plt.subplot(1, 3, 3)
    plt.plot(vyn, vzn, 'm',lw=0.3)
    plt.xlabel('y'); plt.ylabel('z')

    plt.tight_layout()
    plt.savefig(pdf2d, dpi=400, bbox_inches='tight')
    plt.show()
    plt.clf()

    # Atractor completo en 3D1
    pdf3d = pdf2d.replace('.pdf','_3D.pdf')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(vxn, vyn, vzn, 'm', lw=0.2)
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    plt.savefig(pdf3d, dpi=400, bbox_inches='tight')
    plt.show()
    plt.clf()

if __name__ == '__main__':
    main()
