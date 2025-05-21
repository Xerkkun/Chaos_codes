import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def encontrar_valor_mas_cercano(array, valor):
    idx = (np.abs(array - valor)).argmin()
    return idx


def main():
    # Rutas de los archivos y la hoja de GWO
    archivo_de = '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/progress_de_chen_polaco3.txt'
    archivo_pso = '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/progress_pso_chen_polaco2.txt'
    archivo_gwo = '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Chen.xlsx'
    hoja_gwo = 'Semilla 9'

    # Valores específicos que te sirvieron
    valor_de_util = 2.2959 # Ejemplo de valor útil para DE
    valor_pso_util = 2.3138  # Ejemplo de valor útil para PSO
    valor_gwo_util = 2.2936  # Ejemplo de valor útil para GWO

    # Leer datos de DE
    datos_de = []
    with open(archivo_de, 'r') as file:
        for line in file:
            if line.startswith("=") or "n_gen" in line or "Elapsed time" in line:
                continue
            parts = line.split('|')
            value = float(parts[-1].strip())
            datos_de.append(abs(value))  # Tomar el valor absoluto

    # Leer datos de PSO
    datos_pso = []
    with open(archivo_pso, 'r') as file:
        for line in file:
            if line.startswith("=") or "n_gen" in line or "Elapsed time" in line:
                continue
            parts = line.split('|')
            value = float(parts[-1].strip())
            datos_pso.append(abs(value))  # Tomar el valor absoluto

    # Leer datos de GWO desde Excel y ordenar
    df_gwo = pd.read_excel(archivo_gwo, sheet_name=hoja_gwo, engine='openpyxl')
    datos_gwo = df_gwo['Unnamed: 16'].tolist()
    datos_gwo.sort()  # Ordenar los datos de GWO de menor a mayor

    # Convertir listas a arreglos numpy para manipulación
    array_de = np.array(datos_de)
    array_pso = np.array(datos_pso)
    array_gwo = np.array(datos_gwo)

    # Graficar los datos
    # Asumiendo que todas las series tienen la misma longitud
    x = np.arange(len(array_de))

    plt.plot(x, array_de, '-o', label='DE')
    plt.plot(x, array_pso, '-x', label='PSO', color='red')
    plt.plot(x, array_gwo, '-+', label='GWO', color='green')

    # Marcar el valor específico más cercano que te sirvió
    idx_de = encontrar_valor_mas_cercano(array_de, valor_de_util)
    plt.plot(idx_de, array_de[idx_de], 'D',
             color='cyan', markersize=12, label='DE Útil')

    idx_pso = encontrar_valor_mas_cercano(array_pso, valor_pso_util)
    plt.plot(idx_pso, array_pso[idx_pso], 's',
             color='magenta', markersize=12, label='PSO Útil')

    idx_gwo = encontrar_valor_mas_cercano(array_gwo, valor_gwo_util)
    plt.plot(idx_gwo, array_gwo[idx_gwo], 'p',
             color='yellow', markersize=12, label='GWO Útil')

    plt.xlabel('Generations')
    plt.ylabel('D-KY')
    plt.legend()
    plt.savefig('Convergencia_polaco_Chen_resaltada.pdf',
                dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    main()
