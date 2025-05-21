import argparse


def real_to_fixed_point(value, int_bits, frac_bits):
    """
    Convierte un número real a una representación de punto fijo.

    Args:
    value (float): Valor real a convertir.
    int_bits (int): Número de bits para la parte entera incluyendo el bit de signo.
    frac_bits (int): Número de bits para la parte fraccionaria.

    Returns:
    int: Valor en punto fijo como entero.
    """
    # Calcular el factor de escala para la parte fraccionaria
    scale_factor = 2 ** frac_bits

    # Convertir el valor real al entero más cercano que representa el valor de punto fijo
    fixed_point_value = int(round(value * scale_factor))

    # Máscara para manejar el bit de signo y asegurar que no haya desbordamiento
    total_bits = int_bits + frac_bits
    max_value = 2 ** (total_bits - 1) - 1
    min_value = -2 ** (total_bits - 1)

    # Validar el rango del valor convertido
    if fixed_point_value > max_value or fixed_point_value < min_value:
        raise ValueError(
            f"El valor convertido {fixed_point_value} está fuera del rango permitido para {total_bits} bits.")

    # Convertir a representación de dos complementos si es negativo
    if fixed_point_value < 0:
        fixed_point_value = (1 << total_bits) + fixed_point_value

    return fixed_point_value


def main():
    parser = argparse.ArgumentParser(
        description='Convert a real number to fixed-point representation.')
    parser.add_argument('value', type=float,
                        help='The real number to convert.')
    parser.add_argument(
        'int_bits', type=int, help='Number of bits for the integer part, including the sign bit.')
    parser.add_argument('frac_bits', type=int,
                        help='Number of bits for the fractional part.')

    args = parser.parse_args()

    try:
        fixed_point = real_to_fixed_point(
            args.value, args.int_bits, args.frac_bits)
        print(f"Valor en punto fijo (decimal): {fixed_point}")
        print(f"Valor en punto fijo (hexadecimal): {hex(fixed_point)}")
    except ValueError as e:
        print(e)


if __name__ == "__main__":
    main()
