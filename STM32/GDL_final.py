import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Para gráficas 3D

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

def binomial_coef(alpha, mm, decimal):
    """
    Cálculo de coeficientes binomiales para el método GL,
    truncados a 'decimal' dígitos.
    """
    j = np.arange(1, mm + 1, dtype=dtype)
    factors = 1 - (1 + alpha) / j
    c = np.empty(mm + 1, dtype=dtype)
    c[0] = dtype(1.0)
    c[1:] = np.cumprod(factors)
    c = np.trunc(c * (10 ** decimal)) / (10 ** decimal)
    return c.reshape(-1, 1)

def grunwald_letnikov(f, x0, h, alpha, Lm, t_f, system_name, decimal=10):
    """
    Método de Grünwald-Letnikov para derivadas fraccionarias.
    Ahora incluye system_name para nombrar archivos.
    """
    h_alpha = h ** alpha
    m = int(Lm / h)       # memoria en pasos
    k = int(t_f / h)      # número total de iteraciones
    nu = 1 if k < m else m
    mm = k if k < m else m
    d = x0.size           # dimensiones del sistema

    # Inicializar tiempo y estados
    t = np.linspace(0, t_f, k + 1, dtype=dtype)
    x = np.zeros((k + 1, d), dtype=dtype)
    x[0, :] = x0

    # Calcular coeficientes binomiales
    c = binomial_coef(alpha, mm, decimal)

    # Listas para resultados periódicos
    doc_vars = []
    doc_sum  = []

    for i in range(1, k + 1):
        if i > 1:
            idx = np.arange(1, min(nu, i))
            sum_x = np.sum(c[idx] * x[i - idx, :], axis=0)
        else:
            sum_x = np.zeros(d, dtype=dtype)

        # Actualizar estado con fórmula GL
        x[i, :] = f(x[i - 1, :]) * h_alpha - sum_x

        # Guardar cada 100 iteraciones
        if i % 100 == 0:
            vars_line = f"{t[i]:.3f}\t" + "\t".join(f"{val:.10f}" for val in x[i, :])
            sum_line  = "\t".join(f"{val:.15f}" for val in sum_x)
            doc_vars.append(vars_line)
            doc_sum.append(sum_line)

    # Guardar archivos con nombres dinámicos
    np.savetxt(f"GL_{system_name}_python.rnd", doc_vars, fmt="%s")
    np.savetxt(f"sumatoria_GL_{system_name}_python.rnd", doc_sum, fmt="%s")

    return x, t

def main():
    """Función principal: menú y simulación."""
    # Definir menú de sistemas
    tipos = {'1': 'Lorenz', '2': 'Rossler', '3': 'Chen'}
    print("Seleccione sistema a simular:")
    for key, name in tipos.items(): print(f"  {key}: {name}")
    choice = input("Opción (1/2/3): ")
    system_key = tipos.get(choice)
    if system_key is None:
        print("Selección inválida.")
        return

    system_func = {'Lorenz':lorenz_system,'Rossler':rossler_system,'Chen':chen_system}[system_key]


    # Solicitar parámetros al usuario
    alpha   = dtype(float(input("Ingrese orden fraccionario alpha: ")))
    h       = dtype(float(input("Ingrese paso h: ")))
    Lm      = dtype(float(input("Ingrese longitud de memoria Lm: ")))
    t_f     = dtype(float(input("Ingrese tiempo final t_f: ")))
    x0_str = input("Ingrese condiciones iniciales x0 y0 z0 (separadas por espacios): ")
    x0_arr = np.array([float(val) for val in x0_str.strip().split()], dtype=dtype)
    x0, y0, z0 = x0_arr

    # Ejecutar método GL y obtener resultados
    x, t = grunwald_letnikov(system_func, x0, h, alpha, Lm, t_f, system_key)

    # Título y nombre de archivo para gráficas
    pdf2d = f"{system_key}_atractores_EFORK_python.pdf"
    pdf3d = f"{system_key}_atractores_EFORK_3D_python.pdf"


    # Gráficas 2D
    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1); plt.plot(x[:,0], x[:,1], 'm', lw=0.3); plt.xlabel('x'); plt.ylabel('y')
    plt.subplot(1, 3, 2); plt.plot(x[:,0], x[:,2], 'm', lw=0.3); plt.xlabel('x'); plt.ylabel('z')
    plt.subplot(1, 3, 3); plt.plotx(x[:,1], x[:,2], 'm', lw=0.3); plt.xlabel('y'); plt.ylabel('z')
    plt.tight_layout()
    plt.savefig(pdf2d, dpi=400, bbox_inches='tight')
    print(f"Gráfica 2D guardada en: {pdf2d}")
    plt.show()

    # Gráfica 3D
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x[:,0], x[:,1], x[:,2], 'm', lw=0.2)
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    plt.savefig(pdf3d, dpi=400, bbox_inches='tight')
    print(f"Gráfica 3D guardada en: {pdf3d}")
    plt.show()


if __name__ == '__main__':
    main()
