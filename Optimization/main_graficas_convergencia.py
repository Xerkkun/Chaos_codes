import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    
    # Lista para almacenar los valores de la última columna de cada archivo
    all_files_values_de = []
    all_files_values_pso = []
    all_files_values_gwo = []

    xls = pd.ExcelFile('/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Rossler.xlsx' 
, engine='openpyxl')
    nombres_hojas = xls.sheet_names
    
    arreglo_final = []

    for hoja in nombres_hojas:
        #print(hoja)
        df = pd.read_excel('/home/xerk/Escritorio/Códigos para artículo/Artículo PSO ED GW/Analisis de Convergencia Rossler.xlsx' 
, sheet_name=hoja, engine='openpyxl')
        

        # Extraer los datos de la columna 'Unnamed: 16'
        datos_columna = df['Unnamed: 16'].tolist()
        arreglo_final.append(datos_columna)

    # Transponer el arreglo para que cada hoja contribuya a una columna diferente
    arreglo_transpuesto = list(map(list, zip(*arreglo_final)))

    lista_datos_de = []
    lista_datos_pso = []

    for i in range (10):
        # Nombre del archivo basado en el número de iteración
        filename_de = f'progress_de_rossler_polaco{i}.txt'
        filename_pso = f'progress_pso_rossler_polaco{i}.txt'
        
        # Lista para almacenar los valores de la última columna de este archivo
        last_column_values_de = []
        last_column_values_pso = []
        
        # Abrir el archivo de texto y leerlo línea por línea
        with open(filename_de, 'r') as file:
            for line in file:
                # Ignorar líneas que no contienen datos
                if line.startswith("=") or "n_gen" in line or "Elapsed time" in line:
                    continue
                
                # Dividir la línea en sus componentes individuales
                parts = line.split('|')
                
                # Tomar el último valor, eliminar espacios y convertirlo a un número flotante
                value = float(parts[-1].strip())
                
                # Agregar el valor a la lista
                last_column_values_de.append(value)
                
                
        # Abrir el archivo de texto y leerlo línea por línea
        with open(filename_pso, 'r') as file:
            for line in file:
                # Ignorar líneas que no contienen datos
                if line.startswith("=") or "n_gen" in line or "Elapsed time" in line:
                    continue
                
                # Dividir la línea en sus componentes individuales
                parts = line.split('|')
                
                # Tomar el último valor, eliminar espacios y convertirlo a un número flotante
                value = float(parts[-1].strip())
                
                # Agregar el valor a la lista
                last_column_values_pso.append(value)

         # Agregar los valores de este archivo a la lista principal
        all_files_values_de.append(last_column_values_de)
        all_files_values_pso.append(last_column_values_pso)

    # Convertir la lista de listas en un arreglo numpy
    all_files_array_de = -np.array(all_files_values_de)
    all_files_array_pso = -np.array(all_files_values_pso)
    all_files_array_gwo_trans = np.array(arreglo_transpuesto)
    all_files_array_gwo = all_files_array_gwo_trans.T
    
    #Graficar las corridas con los mejores resultados finales para comparar su convergencia
    index_de = np.argmax(all_files_array_de)
    index_pso = np.argmax(all_files_array_pso)
    index_gwo = np.argmax(all_files_array_gwo)
    
    row_de, col_de = np.unravel_index(index_de, all_files_array_de.shape)
    row_pso, col_pso = np.unravel_index(index_pso, all_files_array_pso.shape)
    row_gwo, col_gwo = np.unravel_index(index_gwo, all_files_array_gwo.shape)
    
    x = np.arange(20) #Cantidad de generaciónes 
    
    max_de=all_files_array_de[row_de,col_de]
    max_pso=all_files_array_pso[row_pso,col_pso]
    max_gwo=all_files_array_gwo[row_gwo,col_gwo]
    
    if (max_pso > 2.25) :
        uniques_pso = np.unique(all_files_array_pso)
        second_max_pso = uniques_pso[-2]
        rows_pso, _ = np.where(all_files_array_pso == second_max_pso)
        row_pso = rows_pso[0]
        max_pso=second_max_pso
        
    if (max_de > 2.25):
        uniques_de = np.unique(all_files_array_de)
        second_max_de = uniques_de[-2]
        rows_de, _ = np.where(all_files_array_de == second_max_de)
        row_de = rows_de[0]
        max_de = second_max_de
        
    print(row_de,max_de)
    print(row_pso,max_pso)
    print(row_gwo,max_gwo)
    
    plt.plot(x, all_files_array_de[row_de,:], '-o', label='DE')
    plt.plot(x, all_files_array_pso[row_pso,:], '-x', label='PSO',color='red')
    plt.plot(x, all_files_array_gwo[row_gwo,:], '-+', label='GWO',color='green')    
    plt.xlabel('Generations')
    plt.ylabel('D-KY')
    plt.legend()
    plt.savefig('Convergencia_polaco_rossler.pdf',dpi=300,bbox_inches='tight')
    plt.show()
if __name__=='__main__':
    main()
    
    
