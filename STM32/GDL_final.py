import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Para gráficas 3D

# Seleccionar tipo de dato de mayor precisión disponible
dtype = np.float128 if hasattr(np, 'float128') else np.float64

def lorenz_frac(x):
    """Ecuaciones del sistema caótico de Lorenz."""
    sigma = dtype(10.0)
    beta  = dtype(8.0) / dtype(3.0)
    rho   = dtype(28.0)
    return np.array([
        sigma * (x[1] - x[0]),
        rho * x[0] - x[1] - x[0] * x[2],
        -beta * x[2] + x[0] * x[1]
    ], dtype=dtype)

def rossler_frac(x):
    """Ecuaciones del sistema caótico de Rössler."""
    a = dtype(0.2)
    b = dtype(0.2)
    c = dtype(5.7)
    return np.array([
        -x[1] - x[2],
        x[0] + a * x[1],
        b + x[2] * (x[0] - c)
    ], dtype=dtype)

def chen_frac(x):
    """Ecuaciones del sistema caótico de Chen."""
    u = dtype(7.5)
    v = dtype(1.0)
    w = dtype(5.0)
    return np.array([
        u * (x[1] - x[0]),
        (w - u) * x[0] - x[0] * x[2] + w * x[1],
        x[0] * x[1] - v * x[2]
    ], dtype=dtype)

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

def plot3d(x, y, z, filename):
    """Graficar proyecciones 2D y atractor 3D."""
    plt.figure(figsize=(12, 4))
    # Proyección xy
    plt.subplot(1, 3, 1)
    plt.plot(x, y, 'm', lw=0.3)
    plt.xlabel('x'); plt.ylabel('y')
    # Proyección xz
    plt.subplot(1, 3, 2)
    plt.plot(x, z, 'm', lw=0.3)
    plt.xlabel('x'); plt.ylabel('z')
    # Proyección yz
    plt.subplot(1, 3, 3)
    plt.plot(y, z, 'm',lw=0.3)
    plt.xlabel('y'); plt.ylabel('z')

    plt.tight_layout()
    plt.savefig(filename, dpi=400, bbox_inches='tight')
    plt.show()
    plt.clf()

    # Atractor completo en 3D1
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, 'm', lw=0.2)
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    plt.savefig(filename.replace('.pdf', '_3D.pdf'),dpi=400,bbox_inches= 'tight')
    plt.show()
    plt.clf()

def main():
    """Función principal: menú y simulación."""
    # Definir menú de sistemas
    tipos = {'1': 'Lorenz', '2': 'Rossler', '3': 'Chen'}
    print("Seleccione sistema a simular:")
    for key, name in tipos.items(): print(f"  {key}: {name}")
    choice = input("Opción (1/2/3): ")
    system_name = tipos.get(choice)
    if system_name is None:
        print("Selección inválida.")
        return

    # Mapear función según selección
    f = {
        'Lorenz': lorenz_frac,
        'Rossler': rossler_frac,
        'Chen': chen_frac
    }[system_name]

    # Solicitar parámetros al usuario
    alpha   = dtype(float(input("Ingrese orden fraccionario alpha: ")))
    h       = dtype(float(input("Ingrese paso h: ")))
    Lm      = dtype(float(input("Ingrese longitud de memoria Lm: ")))
    t_f     = dtype(float(input("Ingrese tiempo final t_f: ")))
    x0_str  = input("Ingrese condiciones iniciales separadas por coma: ")
    x0      = np.array([float(val) for val in x0_str.split(',')], dtype=dtype)

    # Ejecutar método GL y obtener resultados
    x, t = grunwald_letnikov(f, x0, h, alpha, Lm, t_f, system_name)

    # Título y nombre de archivo para gráficas
    title = f"{system_name} GL Python α={alpha}"
    filename = f"{system_name}_atractores_GL_python.pdf"

    # Graficar según dimensión
    if x.shape[1] == 3:
        plot3d(x[:,0], x[:,1], x[:,2], filename)
    else:
        print("Gráfica 4D no implementada en este menú.")

if __name__ == '__main__':
    main()
