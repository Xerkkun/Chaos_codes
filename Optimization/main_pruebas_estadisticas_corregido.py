import numpy as np
from scipy.stats import wilcoxon
from scipy.stats import ttest_ind
import pandas as pd

def main():
    
    # Lista para almacenar los valores de la última columna de cada archivo
    all_files_values_de = []
    all_files_values_pso = []
    #all_files_values_gwo = []

    #xls = pd.ExcelFile('/media/programas/INAOE/Doctorado/Optimización/Lorenz/Analisis de Convergencia.xlsx', engine='openpyxl')
    #nombres_hojas = xls.sheet_names
    
    #arreglo_final = []

    # for hoja in nombres_hojas:
    #     #print(hoja)
    #     df = pd.read_excel('/media/programas/INAOE/Doctorado/Optimización/Lorenz/Analisis de Convergencia.xlsx', sheet_name=hoja, engine='openpyxl')
        
    #     # Extraer los datos de la columna 'Unnamed: 4'
    #     datos_columna = df['Unnamed: 4'].tolist()
    #     arreglo_final.append(datos_columna)

    # # Transponer el arreglo para que cada hoja contribuya a una columna diferente
    # arreglo_transpuesto = list(map(list, zip(*arreglo_final)))

    lista_datos_de = []
    lista_datos_pso = []

    for i in range(10):
        nombres_archivos_de = f"last_gen_de_{i}.dat"
        nombres_archivos_pso = f"last_gen_pso_{i}.dat" 
        df_de = pd.read_csv(nombres_archivos_de, sep='\s+', header=None)
        df_pso = pd.read_csv(nombres_archivos_pso, sep='\s+', header=None)
        # Extraer la última columna
        datos_columna_de = df_de.iloc[:, 3].tolist()  # Extraer la cuarta columna usando iloc
        datos_columna_pso = df_pso.iloc[:, 3].tolist()  # Extraer la cuarta columna usando iloc
        lista_datos_de.append(datos_columna_de)
        lista_datos_pso.append(datos_columna_pso)

    # Convertir la lista de listas en un arreglo numpy
    all_files_array_de = -np.array(lista_datos_de)
    all_files_array_pso = -np.array(lista_datos_pso)
    #all_files_array_gwo_trans = np.array(arreglo_transpuesto)
    #all_files_array_gwo = all_files_array_gwo_trans.T
    
    #print(all_files_array_de)

    #print(all_files_array_gwo)
    
    
    # Realizar el test estadístico
    stat = np.zeros(10)
    p = np.zeros(10)
   
    for j in range(9):
        stat[j], p[j] = ttest_ind(all_files_array_pso[j,:], all_files_array_de[j,:])
        # print('---------------------------------------') 
        #print('Estadístico:', stat[j])
        #print('P-valor:', p[j])
    mean_stat = np.mean(stat)
    mean_p = np.mean(p)

    print(mean_stat,mean_p)
    
    # Interpretar el resultado
    alpha = 0.05
    if mean_p > alpha:
        print('No hay evidencia suficiente para rechazar la hipótesis nula (H0), por lo que se asume que las medias son iguales.')
    else:
        print('Se rechaza la hipótesis nula (H0), lo que indica que las medias son diferentes.')
        
        
    mean_de = np.mean(all_files_array_de)
    mean_pso = np.mean(all_files_array_pso)
    #mean_gwo = np.mean(all_files_array_gwo)
    
    #print(mean_de,mean_pso,mean_gwo)
    print(mean_de,mean_pso)
    
    if mean_de > mean_pso:
        print("DE tiene una media mayor que PSO.")
    elif mean_de < mean_pso:
        print("PSO tiene una media mayor que el DE.")
    else:
        print("Ambos grupos tienen la misma media.")
        
    # if mean_de > mean_gwo:
    #     print("DE tiene una media mayor que GWO.")
    # elif mean_de < mean_gwo:
    #     print("GWO tiene una media mayor que el DE.")
    # else:
    #     print("Ambos grupos tienen la misma media.")
        
    # if mean_pso > mean_gwo:
    #     print("PSO tiene una media mayor que GWO.")
    # elif mean_pso < mean_gwo:
    #     print("GWO tiene una media mayor que el PSO.")
    # else:
    #     print("Ambos grupos tienen la misma media.")


    return mean_stat,mean_p


    
if __name__=='__main__':
    main()
    
    