import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt

# Seleccionar tipo de dato de mayor precisión disponible
dtype = np.float128 if hasattr(np, 'float128') else np.float64

def lorenz_system(x, y, z):
    sigma, rho, beta = 10.0, 28.0, 8.0/3.0
    return np.array([sigma*(y-x), x*(rho-z)-y, x*y-beta*z])

def rossler_system(x, y, z):
    a, b, c = 0.2, 0.2, 5.7
    return np.array([-y-z, x+a*y, b+z*(x-c)])

def chen_system(x, y, z):
    u, v, w = 7.5, 1.0, 5.0
    return np.array([u*(y-x), (w-u)*x-x*z+w*y, x*y-v*z])

def memory_fractional(k, t, arr, vtn, h, alpha, nu):
    start = max(0, k-nu)
    sum_ = 0.0
    gamma_term = gamma(2.0 - alpha)
    for j in range(start, k):
        t0 = vtn[j]
        t1 = vtn[j+1]
        v1 = (t - t0)**(1.0 - alpha)
        v2 = (t - t1)**(1.0 - alpha)
        diff = (arr[j+1] - arr[j])
        sum_ += diff * (v1 - v2)
    return sum_ / (h * gamma_term)

def main():
    """Función principal: menú y simulación."""
    # Definir menú de sistemas
    tipos = {'1':'Lorenz','2':'Rossler','3':'Chen'}
    print("Seleccione sistema a simular:")
    for key, name in tipos.items(): print(f"  {key}: {name}")
    choice = input("Opción (1/2/3): ")
    system_key = tipos.get(choice)
    if not system_key:
        print("Selección inválida.")
        return

    system_func = {'Lorenz':lorenz_system,'Rossler':rossler_system,'Chen':chen_system}[system_key]

    # Solicitar parámetros al usuario
    alpha = dtype(float(input("Ingrese orden fraccionario alpha: ")))
    h     = dtype(float(input("Ingrese paso h: ")))
    Lm    = dtype(float(input("Ingrese longitud de memoria Lm: ")))
    t_f   = dtype(float(input("Ingrese tiempo final t_f: ")))
    x0_str = input("Ingrese condiciones iniciales x0 y0 z0 (separadas por espacios): ")
    x0_arr = np.array([float(val) for val in x0_str.strip().split()], dtype=dtype)
    x0, y0, z0 = x0_arr

    # Parámetros EFORK
    N1 = int(np.ceil(t_f / h))
    nu = int(Lm / h)
    ha = h ** alpha

    gamma1 = gamma(1.0 + alpha)
    gamma2 = gamma(1.0 + 2 * alpha)
    gamma3 = gamma(1.0 + 3 * alpha)
    c2 = (1.0 / (2.0 * gamma1)) ** (1.0 / alpha)
    c3 = (1.0 / (4.0 * gamma1)) ** (1.0 / alpha)
    a21 = 1.0 / (2.0 * gamma1 * gamma1)
    a31 = (gamma1**2 * gamma(2*alpha+1) + 2.0*gamma(2*alpha+1)**2 - gamma(3*alpha+1)) / \
          (4.0 * gamma1**2 * (2.0 * gamma(2*alpha+1)**2 - gamma(3*alpha+1)))
    a32 = - gamma(2*alpha+1) / (4.0 * (2.0 * gamma(2*alpha+1)**2 - gamma(3*alpha+1)))
    w1 = (8.0 * gamma1**3 * gamma(1.0 + 2 * alpha)**2 - 6.0 * gamma1**3 * gamma(1.0 + 3 * alpha) + gamma(1.0 + 2 * alpha) * gamma(1.0 + 3 * alpha)) / \
         (gamma1 * gamma(1.0 + 2 * alpha) * gamma(1.0 + 3 * alpha))
    w2 = 2.0 * gamma1**2 * (4.0 * gamma(1.0 + 2 * alpha)**2 - gamma(1.0 + 3 * alpha)) / \
         (gamma(1.0 + 2 * alpha) * gamma(1.0 + 3 * alpha))
    w3 = -8.0 * gamma1**2 * (2.0 * gamma(1.0 + 2 * alpha)**2 - gamma(1.0 + 3 * alpha)) / \
         (gamma(1.0 + 2 * alpha) * gamma(1.0 + 3 * alpha))

    # Inicializar arrays
    vtn = np.zeros(N1+1)
    vxn = np.zeros(N1+1)
    vyn = np.zeros(N1+1)
    vzn = np.zeros(N1+1)
    vtn[0], vxn[0], vyn[0], vzn[0] = 0.0, x0, y0, z0

    # Primer paso (n=0): SIN MEMORIA
    x_n, y_n, z_n = x0, y0, z0
    dx, dy, dz = system_func(x_n, y_n, z_n)
    K1x, K1y, K1z = ha * dx, ha * dy, ha * dz
    dx2, dy2, dz2 = system_func(x_n + a21*K1x, y_n + a21*K1y, z_n + a21*K1z)
    K2x, K2y, K2z = ha * dx2, ha * dy2, ha * dz2
    dx3, dy3, dz3 = system_func(x_n + a31*K2x + a32*K1x,
                                y_n + a31*K2y + a32*K1y,
                                z_n + a31*K2z + a32*K1z)
    K3x, K3y, K3z = ha * dx3, ha * dy3, ha * dz3

    x_n1 = x_n + w1*K1x + w2*K2x + w3*K3x
    y_n1 = y_n + w1*K1y + w2*K2y + w3*K3y
    z_n1 = z_n + w1*K1z + w2*K2z + w3*K3z

    vtn[1], vxn[1], vyn[1], vzn[1] = h, x_n1, y_n1, z_n1
    x_n, y_n, z_n = x_n1, y_n1, z_n1

    # Ciclo principal (n>0)
    for n in range(1, N1):
        tn = n * h
        mem_x = memory_fractional(n, tn, vxn, vtn, h, alpha, nu)
        mem_y = memory_fractional(n, tn, vyn, vtn, h, alpha, nu)
        mem_z = memory_fractional(n, tn, vzn, vtn, h, alpha, nu)

        dx, dy, dz = system_func(x_n, y_n, z_n)
        f1x, f1y, f1z = dx - mem_x, dy - mem_y, dz - mem_z
        K1x, K1y, K1z = ha * f1x, ha * f1y, ha * f1z

        dx2, dy2, dz2 = system_func(x_n + a21*K1x, y_n + a21*K1y, z_n + a21*K1z)
        K2x, K2y, K2z = ha * dx2, ha * dy2, ha * dz2

        dx3, dy3, dz3 = system_func(x_n + a31*K2x + a32*K1x,
                                    y_n + a31*K2y + a32*K1y,
                                    z_n + a31*K2z + a32*K1z)
        K3x, K3y, K3z = ha * dx3, ha * dy3, ha * dz3

        x_n1 = x_n + w1*K1x + w2*K2x + w3*K3x
        y_n1 = y_n + w1*K1y + w2*K2y + w3*K3y
        z_n1 = z_n + w1*K1z + w2*K2z + w3*K3z

        vtn[n+1], vxn[n+1], vyn[n+1], vzn[n+1] = (n+1)*h, x_n1, y_n1, z_n1
        x_n, y_n, z_n = x_n1, y_n1, z_n1

    # Guardar resultados
    datos = np.column_stack((vtn, vxn, vyn, vzn))
    out_txt = f"EFORK_{system_key}_c.txt"
    np.savetxt(out_txt, datos, delimiter='\t', header='t\tx\ty\tz', comments='')
    print(f"Resultados guardados en: {out_txt}")

    # Nombres de PDF
    pdf2d = f"{system_key}_atractores_EFORK_python.pdf"
    pdf3d = f"{system_key}_atractores_EFORK_3D_python.pdf"

    # Gráficas 2D
    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1); plt.plot(vxn, vyn, 'm', lw=0.3); plt.xlabel('x'); plt.ylabel('y')
    plt.subplot(1, 3, 2); plt.plot(vxn, vzn, 'm', lw=0.3); plt.xlabel('x'); plt.ylabel('z')
    plt.subplot(1, 3, 3); plt.plot(vyn, vzn, 'm', lw=0.3); plt.xlabel('y'); plt.ylabel('z')
    plt.tight_layout()
    plt.savefig(pdf2d, dpi=400, bbox_inches='tight')
    print(f"Gráfica 2D guardada en: {pdf2d}")
    plt.show()

    # Gráfica 3D
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(vxn, vyn, vzn, 'm', lw=0.2)
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    plt.savefig(pdf3d, dpi=400, bbox_inches='tight')
    print(f"Gráfica 3D guardada en: {pdf3d}")
    plt.show()

if __name__ == '__main__':
    main()
