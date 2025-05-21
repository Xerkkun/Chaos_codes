import numpy as np
import matplotlib.pyplot as plt


class FixedPoint:
    def __init__(self, value, total_bits, integer_bits):
        self.total_bits = total_bits
        self.integer_bits = integer_bits
        self.fractional_bits = total_bits - integer_bits - 1
        self.max_val = (1 << (total_bits - 1)) - 1
        self.min_val = -(1 << (total_bits - 1))
        self.scale = 1 << self.fractional_bits
        self.value = self.float_to_fixed(value)

    def float_to_fixed(self, value):
        fixed_value = int(round(value * self.scale))
        if fixed_value > self.max_val:
            fixed_value = self.max_val
        elif fixed_value < self.min_val:
            fixed_value = self.min_val
        return fixed_value

    def fixed_to_float(self, fixed_value):
        return fixed_value / self.scale

    def __add__(self, other):
        return FixedPoint(self.fixed_to_float(self.value + other.value), self.total_bits, self.integer_bits)

    def __sub__(self, other):
        return FixedPoint(self.fixed_to_float(self.value - other.value), self.total_bits, self.integer_bits)

    def __mul__(self, other):
        result = (self.value * other.value) >> self.fractional_bits
        return FixedPoint(self.fixed_to_float(result), self.total_bits, self.integer_bits)

    def __truediv__(self, other):
        result = (self.value << self.fractional_bits) // other.value
        return FixedPoint(self.fixed_to_float(result), self.total_bits, self.integer_bits)

    def to_float(self):
        return self.fixed_to_float(self.value)

# Define las ecuaciones del sistema de Lorenz en punto fijo


def lorenz_equations_fixed(X, sigma, rho, beta):
    x, y, z = X
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return dx, dy, dz

# Método de Euler hacia adelante con punto fijo


def euler_forward_fixed(X, h, sigma, rho, beta):
    dx, dy, dz = lorenz_equations_fixed(X, sigma, rho, beta)
    X = [X[0] + h * dx, X[1] + h * dy, X[2] + h * dz]
    return X

# Función principal para simular el sistema de Lorenz


def simulate_lorenz_fixed(h, iterations, transient, sigma, rho, beta, total_bits=32, integer_bits=13):
    # Convertir los parámetros a punto fijo
    sigma = FixedPoint(sigma, total_bits, integer_bits)
    rho = FixedPoint(rho, total_bits, integer_bits)
    beta = FixedPoint(beta, total_bits, integer_bits)
    h = FixedPoint(h, total_bits, integer_bits)

    # Condiciones iniciales
    x0 = FixedPoint(-10.261005765339908, total_bits, integer_bits)
    y0 = FixedPoint(10.261005765339908, total_bits, integer_bits)
    z0 = FixedPoint(134.280424986204, total_bits, integer_bits)
    X = [x0, y0, z0]

    # Simular el sistema
    sol = np.zeros((iterations + 1, 3))
    sol[0] = [x0.to_float(), y0.to_float(), z0.to_float()]

    for i in range(iterations):
        X = euler_forward_fixed(X, h, sigma, rho, beta)
        sol[i + 1] = [X[0].to_float(), X[1].to_float(), X[2].to_float()]

        # Condición de paro por desbordamiento
        if abs(X[0].to_float()) > 1E3 or abs(X[1].to_float()) > 1E3 or abs(X[2].to_float()) > 1E3:
            break

    # Quitar el transitorio
    sol = sol[transient:, :]

    # Graficar los atractores de Lorenz
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(sol[:, 0], sol[:, 1], sol[:, 2], lw=0.5)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()


# Parámetros de entrada
h = 0.00894
iterations = 1000000
transient = 50000
sigma, rho, beta = 5.301300510626442808e+00, 1.352804249862039967e+02, 7.840922407503264635e-01

# Ejecutar la simulación
simulate_lorenz_fixed(h, iterations, transient, sigma, rho, beta)
