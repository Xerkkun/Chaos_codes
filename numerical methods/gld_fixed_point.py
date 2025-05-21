import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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

    def to_float(self):
        return self.fixed_to_float(self.value)

    def __neg__(self):
        return FixedPoint(-self.to_float(), self.total_bits, self.integer_bits)

    def __add__(self, other):
        if isinstance(other, FixedPoint):
            return FixedPoint(self.to_float() + other.to_float(), self.total_bits, self.integer_bits)
        else:
            return FixedPoint(self.to_float() + other, self.total_bits, self.integer_bits)

    def __sub__(self, other):
        if isinstance(other, FixedPoint):
            return FixedPoint(self.to_float() - other.to_float(), self.total_bits, self.integer_bits)
        else:
            return FixedPoint(self.to_float() - other, self.total_bits, self.integer_bits)

    def __mul__(self, other):
        if isinstance(other, FixedPoint):
            result = (self.value * other.value) >> self.fractional_bits
            return FixedPoint(self.fixed_to_float(result), self.total_bits, self.integer_bits)
        else:
            result = (self.value * other) >> self.fractional_bits
            return FixedPoint(self.fixed_to_float(result), self.total_bits, self.integer_bits)

    def __truediv__(self, other):
        if isinstance(other, FixedPoint):
            result = (self.value << self.fractional_bits) // other.value
            return FixedPoint(self.fixed_to_float(result), self.total_bits, self.integer_bits)
        else:
            result = (self.value << self.fractional_bits) // other
            return FixedPoint(self.fixed_to_float(result), self.total_bits, self.integer_bits)

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        return FixedPoint(other, self.total_bits, self.integer_bits).__sub__(self)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __rtruediv__(self, other):
        return FixedPoint(other, self.total_bits, self.integer_bits).__truediv__(self)


def lorenz_frac_fixed(x, sigma, beta, rho):
    xdot = np.zeros_like(x, dtype=object)
    xdot[0] = sigma * (x[1] - x[0])
    xdot[1] = (rho * x[0]) - x[1] - (x[0] * x[2])
    xdot[2] = (-beta * x[2]) + (x[0] * x[1])
    return xdot


def binomial_coef_fixed(alpha, mm, total_bits, integer_bits):
    c = [FixedPoint(1.0, total_bits, integer_bits)]
    with open("binomial_coefficients.txt", "w") as file:
        file.write(
            f"{c[0].to_float()} (decimal), {c[0].value} (fixed point)\n")
        for j in range(1, mm + 1):
            prev_coef = c[j - 1]
            new_coef = FixedPoint(prev_coef.to_float() *
                                  (1 + alpha) / j, total_bits, integer_bits)
            c.append(prev_coef - new_coef)
            file.write(
                f"{c[j].to_float()} (decimal), {c[j].value} (fixed point)\n")
    return c


def grunwald_letnikov_fixed(x, h, h_alpha, k, alpha, x_t, nu, d, mm, sigma, beta, rho, total_bits, integer_bits):
    sum_x = np.array([FixedPoint(0.0, total_bits, integer_bits)
                     for _ in range(d)])
    c = binomial_coef_fixed(alpha, mm, total_bits, integer_bits)
    new_x = np.zeros(d, dtype=object)
    aux_x = np.zeros((nu, d), dtype=object)
    aux_x[0, :] = x[0, :]

    for i in range(1, k + 1):
        if nu == 1:
            for j in range(1, i + 1):
                for idx in range(d):
                    sum_x[idx] = sum_x[idx] + (c[j] * x[i - j, idx])
            x[i, :] = lorenz_frac_fixed(
                x[i - 1, :], sigma, beta, rho) * h_alpha - sum_x
        else:
            for j in range(1, nu + 1):
                for idx in range(d):
                    sum_x[idx] = sum_x[idx] + (c[j] * aux_x[nu - j, idx])
            if i < mm:
                x[i, :] = lorenz_frac_fixed(
                    x[i - 1, :], sigma, beta, rho) * h_alpha - sum_x
                aux_x[i] = x[i, :]
            else:
                new_x = lorenz_frac_fixed(
                    x[i - 1, :], sigma, beta, rho) * h_alpha - sum_x
                aux_x = np.concatenate((aux_x[1:], [new_x]), axis=0)
                x[i, :] = new_x

        sum_x = np.array([FixedPoint(0.0, total_bits, integer_bits)
                         for _ in range(d)])
    return x


def main():
    total_bits = 32
    integer_bits = 10

    alpha = 0.98
    sigma = FixedPoint(10.0, total_bits, integer_bits)
    beta = FixedPoint(8. / 3., total_bits, integer_bits)
    rho = FixedPoint(28.0, total_bits, integer_bits)
    x_0 = np.array([FixedPoint(0.1, total_bits, integer_bits),
                    FixedPoint(0.1, total_bits, integer_bits),
                    FixedPoint(0.1, total_bits, integer_bits)], dtype=object)

    t_0 = 0.0
    t_f = 1000.0
    h = FixedPoint(0.01, total_bits, integer_bits)
    h_alpha = FixedPoint(h.to_float() ** alpha, total_bits, integer_bits)

    Lm = 0.1
    m = int(Lm / h.to_float())
    k = int((t_f - t_0) / h.to_float())

    if k < m:
        nu, mm = 1, k
    else:
        nu, mm = m, m

    d = 3
    x = np.zeros((k + 1, d), dtype=object)
    t = np.linspace(t_0, t_f, len(x))
    x_t = np.zeros((k + 1, 1))
    x_t[:, 0] = t
    x[0, :] = x_0

    x = grunwald_letnikov_fixed(
        x, h, h_alpha, k, alpha, x_t, nu, d, mm, sigma, beta, rho, total_bits, integer_bits)

    t = x_t[:, 0]
    xx = np.array([xi.to_float() for xi in x[:, 0]])
    y = np.array([yi.to_float() for yi in x[:, 1]])
    z = np.array([zi.to_float() for zi in x[:, 2]])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xx, y, z, lw=0.5)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()


if __name__ == '__main__':
    main()
