import math
from scipy.special import gamma
# ================================================================================================================
def binomial_derivative_gamma(n, k):
    """Calcula el coeficiente binomial para la derivada fraccionaria de orden n"""
    return gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1))

def binomial_integral_gamma(n, k):
    """Calcula el coeficiente binomial para la integral fraccionaria de orden n"""
    return (-1)**k * gamma(n + k) / (gamma(k + 1) * gamma(n))

def calculate_coefficients(order, num_coefficients):
    """Calcula los coeficientes para la derivada e integral fraccionaria"""
    derivative_coeffs = [binomial_derivative_gamma(order, k) for k in range(num_coefficients)]
    integral_coeffs = [binomial_integral_gamma(order, k) for k in range(num_coefficients)]
    return derivative_coeffs, integral_coeffs
#================================================================================================================
def calculate_binomial_coefficients(alpha, num_coefficients):
    """Calcula los coeficientes binomiales fraccionarios utilizando la expresión iterativa."""
    coefficients = [1.0]  # C_alpha_0 = 1

    for j in range(1, num_coefficients):
        prev_coeff = coefficients[j-1]
        new_coeff = (1 - (1 + alpha) / j) * prev_coeff
        coefficients.append(new_coeff)
    return coefficients

def calculate_binomial_integral_coefficients(alpha, num_coefficients):
    """Calcula los coeficientes binomiales fraccionarios para la integral utilizando la expresión iterativa."""
    coefficients = [1.0]  # C_alpha_0 = 1

    for j in range(1, num_coefficients):
        prev_coeff = coefficients[j-1]
        # new_coeff = (-1)**j * (1 - (1 - alpha) / j) * prev_coeff
        new_coeff = (1 - (1 - alpha) / j) * prev_coeff
        coefficients.append(new_coeff)
    return coefficients
# ================================================================================================================
# Ejemplo de uso:
alpha = 0.5
num_coefficients = 10
coefficients = calculate_binomial_integral_coefficients(
    alpha, num_coefficients)

print(
    f"Coeficientes binomiales fraccionarios para la integral con alpha = {alpha}:")
for j, coeff in enumerate(coefficients):
    print(f"  C_alpha_{j} = {coeff:.6f}")

# Ejemplo de uso:
alpha = 0.5
num_coefficients = 10
coefficients = calculate_binomial_coefficients(alpha, num_coefficients)

print(f"Coeficientes binomiales fraccionarios para alpha = {alpha}:")
for j, coeff in enumerate(coefficients):
    print(f"  C_alpha_{j} = {coeff:.6f}")

# Ejemplo de uso:
order = 0.5
num_coefficients = 10
derivative_coeffs, integral_coeffs = calculate_coefficients(order, num_coefficients)

print(f"Coeficientes de la derivada fraccionaria de orden {order}:")
for i, coeff in enumerate(derivative_coeffs):
    print(f"  C({order}, {i}) = {coeff:.6f}")

print(f"\nCoeficientes de la integral fraccionaria de orden {order}:")
for i, coeff in enumerate(integral_coeffs):
    print(f"  C(-{order}, {i}) = {coeff:.6f}")
