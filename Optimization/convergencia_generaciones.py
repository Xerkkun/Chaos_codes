import matplotlib.pyplot as plt


def leer_archivo_dat_y_extraer_f_min_positivos(ruta_archivo):
    valores_f_min_positivos = []
    with open(ruta_archivo, 'r') as archivo:
        # Omitir encabezado y pie de página
        for linea in archivo.readlines()[2:-1]:
            partes = linea.split('|')
            if len(partes) > 3:  # Asegurar que la línea tiene suficientes columnas
                f_min = partes[-1].strip()  # Obtener y limpiar el último valor
                try:
                    # Convertir a valor absoluto y agregar a la lista
                    valores_f_min_positivos.append(abs(float(f_min)))
                except ValueError:
                    pass  # Manejar posibles errores de conversión
    return valores_f_min_positivos


def graficar_comparacion_positivos(archivo1, archivo2):
    valores1 = leer_archivo_dat_y_extraer_f_min_positivos(archivo1)
    valores2 = leer_archivo_dat_y_extraer_f_min_positivos(archivo2)

    # Crear una figura y un eje para la gráfica
    plt.figure(figsize=(10, 6))
    plt.plot(valores1, label='Archivo 1', marker='o')
    plt.plot(valores2, label='Archivo 2', marker='x')

    # Añadir título y etiquetas a los ejes
    plt.title('Comparación de f_min positivos entre dos archivos')
    plt.xlabel('Número de línea')
    plt.ylabel('f_min (valores positivos)')

    # Añadir una leyenda
    plt.legend()

    # Mostrar la gráfica
    plt.show()


# Ejemplo de uso
archivo1 = '/home/xerk/Escritorio/Códigos para artículo/progress_de_polaco0.txt'
archivo2 = '/home/xerk/Escritorio/Códigos para artículo/progress_de_polaco9.txt'

graficar_comparacion_positivos(archivo1, archivo2)
