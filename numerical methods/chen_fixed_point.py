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

    def __neg__(self):
        return FixedPoint(self.fixed_to_float(-self.value), self.total_bits, self.integer_bits)

    def to_float(self):
        return self.fixed_to_float(self.value)

# Define las ecuaciones del sistema de Chen en punto fijo


def chen_equations_fixed(X, a, b, c):
    x, y, z = X
    dx = a * (y - x)
    dy = (c - a) * x - x * z + c * y
    dz = x * y - b * z
    return dx, dy, dz

# Método de Euler hacia adelante con punto fijo


def euler_forward_fixed(X, h, a, b, c):
    dx, dy, dz = chen_equations_fixed(X, a, b, c)
    X = [X[0] + h * dx, X[1] + h * dy, X[2] + h * dz]
    return X

# Función principal para simular el sistema de Chen


def simulate_chen_fixed(h, iterations, transient, a, b, c, total_bits=32, integer_bits=12):
    # Convertir los parámetros a punto fijo
    a = FixedPoint(a, total_bits, integer_bits)
    b = FixedPoint(b, total_bits, integer_bits)
    c = FixedPoint(c, total_bits, integer_bits)
    h = FixedPoint(h, total_bits, integer_bits)

    # Condiciones iniciales
    x0 = FixedPoint(4.246, total_bits, integer_bits)
    y0 = FixedPoint(4.728, total_bits, integer_bits)
    z0 = FixedPoint(13.470, total_bits, integer_bits)
    X = [x0, y0, z0]

    # Simular el sistema
    sol = np.zeros((iterations + 1, 3))
    sol[0] = [x0.to_float(), y0.to_float(), z0.to_float()]

    for i in range(iterations):
        X = euler_forward_fixed(X, h, a, b, c)
        sol[i + 1] = [X[0].to_float(), X[1].to_float(), X[2].to_float()]

        # Condición de paro por desbordamiento
        if abs(X[0].to_float()) > 1E3 or abs(X[1].to_float()) > 1E3 or abs(X[2].to_float()) > 1E3:
            break

    # Quitar el transitorio
    sol = sol[transient:, :]

    # Graficar los atractores de Chen
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(sol[:, 0], sol[:, 1], sol[:, 2], lw=0.5)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()


# Parámetros de entrada
h = 0.00193
iterations = 100000
transient = 50000
a, b, c = 2.821888891274968714e+01,  7.066636808327453334e-01,  2.475217691457946856e+01



# Ejecutar la simulación
simulate_chen_fixed(h, iterations, transient, a, b, c)
