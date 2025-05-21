import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    # Rutas de los archivos y la hoja de GWO
    archivo_de = '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/progress_de_chen_polaco3.txt'
    archivo_pso = '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/progress_pso_chen_polaco2.txt'
    archivo_gwo = '/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Chen.xlsx'
    hoja_gwo = 'Semilla 9'

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
    plt.xlabel('Generations')
    plt.ylabel('D-KY')
    plt.legend()
    plt.savefig('Convergencia_polaco_Chen.pdf', dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    main()
