"""
   Información
   Programa que utiliza la definición del operador diferencial de Grünwald-Letnikov
   para aproximar la solución numérica de un sistema de ecuaciones de orden
   fraccionario alpha con valores iniciales.
   Se utiliza el sistema caótico de Lorenz para calibrar el método numérico.
"""
import numpy as np
import matplotlib.pyplot as plt

# Seleccionar el tipo de dato de mayor precisión disponible
dtype = np.float128 if hasattr(np, 'float128') else np.float64

def lorenz_frac(x):
    """Ecuaciones del sistema caótico de Lorenz, con alta precisión."""

    # alpha=0.98
    # h=0.01
    # Lm=1

    sigma = dtype(10.0)
    beta = dtype(8.0) / dtype(3.0)
    rho = dtype(28.0)
    return np.array([
        sigma * (x[1] - x[0]),
        rho * x[0] - x[1] - x[0] * x[2],
        -beta * x[2] + x[0] * x[1]
    ], dtype=dtype)

def rossler_frac(x):
    """Ecuaciones del sistema caótico de Rössler, con alta precisión."""

    # alpha=0.985
    # h=0.01
    # Lm=1

    a = dtype(0.2)
    b = dtype(0.2)
    c = dtype(5.7)
    return np.array([
        - x[1] - x[2],
        x[0] + a * x[1],
        b + x[2]*(x[0] - c)
    ], dtype=dtype)

def chen_frac(x):
    """Ecuaciones del sistema caótico de Chen, con alta precisión."""
    # alpha=0.91
    # h=0.01
    # Lm=1

    u = dtype(7.5)
    v = dtype(1.0)
    w = dtype(5.0)
    return np.array([
        u*(x[1] - x[0]),
        (w - u)*x[0] - x[0]*x[2] + w*x[1],
        x[0]*x[1] - v*x[2]
    ], dtype=dtype)

def ho2_system(x):
    """Sistema HO2 (no se utiliza en este ejemplo)."""
    b = dtype(0.1)
    a = dtype(0.2)
    xdot = np.zeros_like(x, dtype=dtype)
    xdot[0] = x[1] * x[2]
    xdot[1] = x[0] - x[1] - a * x[3]
    xdot[2] = 1 - x[0] ** 2
    xdot[3] = b * x[1]
    return xdot

def binomial_coef(alpha, mm, decimal):
    """
    Cálculo vectorizado de los coeficientes binomiales con truncamiento decimal.
    Se computa:
      c[0] = 1
      c[j] = producto_{i=1}^{j} (1 - (1+alpha)/i)
    """
    j = np.arange(1, mm + 1, dtype=dtype)
    factors = 1 - (1 + alpha) / j
    c = np.empty(mm + 1, dtype=dtype)
    c[0] = dtype(1.0)
    c[1:] = np.cumprod(factors)
    factor = 10 ** decimal
    c = np.trunc(c * factor) / factor
    with open("lorenz_coeficientes.rnd", "w") as arch:
        for val in c:
            arch.write(f"{val:.15f}\n")
    return c.reshape(-1, 1)

def grunwald_letnikov(x, h, h_alpha, k, alpha, x_t, nu, d, mm, decimal, m):
    """
    Método de Grünwald-Letnikov para la derivada fraccionaria.
    Se vectoriza la suma interna y se bufferiza la escritura a disco para mejorar el rendimiento.
    """
    c = binomial_coef(alpha, mm, decimal)
    variables_lines = []
    sumatoria_lines = []
    
    for i in range(1, k + 1):
        # Calcular la suma de coeficientes multiplicados por los estados anteriores
        if i > 1:
            idx = np.arange(1, min(nu, i) if nu != 1 else i)
            sum_x = np.sum(c[idx] * x[i - idx, :], axis=0)
        else:
            sum_x = np.zeros(d, dtype=dtype)
            
        # Se asigna el nuevo estado
        x[i, :] = lorenz_frac(x[i - 1, :]) * h_alpha - sum_x
        
        # Almacenar datos cada 100 iteraciones
        if d == 3:
            if i % 100 == 0:
                print(f"{x_t[i,0]:.3f} {x[i,0]:.10f} {x[i,1]:.10f} {x[i,2]:.10f}")
            variables_lines.append(f"{x_t[i,0]:.3f}\t{x[i,0]:.10f}\t{x[i,1]:.10f}\t{x[i,2]:.10f}")
            sumatoria_lines.append(f"{sum_x[0]:.15f}\t{sum_x[1]:.15f}\t{sum_x[2]:.15f}")
        else:
            if i % 100 == 0:
                print(f"{x_t[i,0]:.3f} {x[i,0]:.10f} {x[i,1]:.10f} {x[i,2]:.10f} {x[i,3]:.10f}")
            variables_lines.append(f"{x_t[i,0]:.3f}\t{x[i,0]:.10f}\t{x[i,1]:.10f}\t{x[i,2]:.10f}\t{x[i,3]:.10f}")
            sumatoria_lines.append(f"{sum_x[0]:.15f}\t{sum_x[1]:.15f}\t{sum_x[2]:.15f}\t{sum_x[3]:.15f}")
    
    with open("lorenz_variables.rnd", "w") as arch1:
        arch1.write("\n".join(variables_lines))
    with open("lorenz_sumatoria.rnd", "w") as arch2:
        arch2.write("\n".join(sumatoria_lines))
    
    return x

def grafica3d(x, y, z, t):
    """Graficar atractores en 3D, guardar en PDF y mostrar la gráfica."""
    plt.figure (figsize=(12, 4))
    plt.subplot(1, 3, 1)
    plt.plot(x, y, "m", lw=0.3)
    plt.xlabel("x")
    plt.ylabel("y")
    
    plt.subplot(1, 3, 2)
    plt.plot(x, z, "m", lw=0.3)
    plt.xlabel("x")
    plt.ylabel("z")
    
    plt.subplot(1, 3, 3)
    plt.plot(y, z, "m", lw=0.3)
    plt.xlabel("y")
    plt.ylabel("z")
    
    plt.tight_layout()
    plt.savefig("lorenz_atractores.pdf", dpi=400, bbox_inches='tight')
    plt.show()  # Se muestra la gráfica al final
    plt.clf()

def grafica4d(x, y, z, w, t):
    """Graficar atractores en 4D, guardar en PDF y mostrar la gráfica."""
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(x, y, "m", lw=0.3)
    plt.xlabel("x")
    plt.ylabel("y")
    
    plt.subplot(2, 2, 2)
    plt.plot(y, z, "m", lw=0.3)
    plt.xlabel("y")
    plt.ylabel("z")
    
    plt.subplot(2, 2, 3)
    plt.plot(x, z, "m", lw=0.3)
    plt.xlabel("x")
    plt.ylabel("z")
    
    plt.subplot(2, 2, 4)
    plt.plot(x, w, "m", lw=0.3)
    plt.xlabel("x")
    plt.ylabel("w")
    
    plt.tight_layout()
    plt.savefig("lorenz_atractores_2.pdf", dpi=300, bbox_inches='tight')
    plt.show()  # Se muestra la gráfica al final
    plt.clf()

def main():
    """Integración numérica del sistema de Lorenz con el método de Grünwald-Letnikov,
    optimizado en tiempo y precisión."""
    decimal = 10  # Precisión decimal deseada
    alpha = dtype(0.98)  # Orden fraccionario
    x_0 = np.array([0.1, 0.1, 0.1], dtype=dtype)  # Condición inicial
    
    t_0 = dtype(0.0)
    t_f = dtype(250.0)
    h = dtype(0.01)
    h_alpha = h ** alpha
    
    Lm = dtype(1)         # Longitud de memoria
    m = int(Lm / h)         # Número máximo de iteraciones en memoria
    k = int((t_f - t_0) / h)  # Número total de puntos (por ejemplo, 10000)
    
    if k < m:
        nu = 1
        mm = k
    else:
        nu = m
        mm = m
    
    d = 3  # Dimensiones del sistema (3 para Lorenz, 4 para HO2)
    x = np.zeros((k + 1, d), dtype=dtype)
    t = np.linspace(t_0, t_f, k + 1).astype(dtype)
    x_t = t.reshape(-1, 1)
    
    x[0, :] = x_0
    x = grunwald_letnikov(x, h, h_alpha, k, alpha, x_t, nu, d, mm, decimal, m)
    
    if d == 3:
        xx, y, z = x[:, 0], x[:, 1], x[:, 2]
        grafica3d(xx, y, z, t)
    else:
        xx, y, z, w = x[:, 0], x[:, 1], x[:, 2], x[:, 3]
        grafica4d(xx, y, z, w, t)

if __name__ == '__main__':
    main()
