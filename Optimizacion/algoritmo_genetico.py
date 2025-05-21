import random
import numpy as np
import math
import matplotlib.pyplot as plt
#===============================================================================
#Funciones
#===============================================================================

#def generar_poblacion_inicial(tamano_poblacion, longitud_cromosoma):
#    """se generan individuos de la población de tamaño tamano_poblacion 
#    con una longitud de cromosoma de longitud_cromosoma. Cada cromosoma 
#    es una lista de valores binarios generados aleatoriamente."""
#    poblacion = []
#    for _ in range(tamano_poblacion):
#        cromosoma = [random.randint(0, 1) for _ in range(longitud_cromosoma)]
#        poblacion.append(cromosoma)
#    return poblacion

def fx(x):
    return -(0.1+(1-x)**2-0.1*math.cos(6*math.pi*(1-x)))+2

def listToDecimal(num):
    decimal=0
    for i in range(len(num)):
        decimal+=num[i]*10**(-i)
    return decimal

def mutate(individuals, prob, pool):
    for i in range(len(individuals)):
        mutate_individual=individuals[i]
        if np.random.random() < prob:
            mutation = np.random.choice(pool[0])
            mutate_individual = [mutation] + mutate_individual[1:]
        
        for j in range(1,len(mutate_individual)):
            if np.random.random() < prob:
                mutation = np.random.choice(pool[1])
                mutate_individual = mutate_individual[0:j] + [mutation] + mutate_individual[j+1:]
        individuals[i] = mutate_individual

def main():
    y_axis = []
    x_axis=np.arange(0,2,0.02)

    for num in x_axis:
        y_axis.append(fx(num))
    
    #plt.plot(x_axis,y_axis)
    #plt.show()

    # Parámetros del algoritmo genético
    #tamaño_poblacion = 100
    #probabilidad_cruce = 0.8
    #probabilidad_mutacion = 0.1
    #numero_generaciones = 10
    #longitud_cromosoma = 8
    
    #poblacion = generar_poblacion_inicial(tamaño_poblacion, longitud_cromosoma)
    #print(str(poblacion) + '\n')
    
    #individuo de ejemplo 
    #x = 0.54
    #y = fx(x)
    
    #plt.plot(x_axis,y_axis)
    #plt.plot(x,y,'x')
    #plt.show()

    x = [0,5,4]
    
    # Necesitamos representar los genes como una lista para poder realizar la 
    # mutación y el entrecruzamiento. Así que debemos manejar una función 
    # que convierta una lista en números decimales.
    
    #decimal = listToDecimal(x)
    
    # Me gustaría manejar individuos con un código genético más grande, 
    # así que voy a definir un ind_size con el cual generaré un primer 
    # código genético aleatorio. 
    
    ind_size = 15
    #Genetic pool
    genetic_pool=[[0,1],[0,1,2,3,4,5,6,7,8,9]]
    
    # Se crea una lista vacía individuo. Luego, se agrega un elemento aleatorio 
    # de genetic_pool[0] utilizando np.random.choice y el método append. 
    # A continuación, se agregan elementos adicionales a individuo utilizando 
    # np.random.choice con genetic_pool[1] y ind_size - 1 repeticiones, y se 
    # concatenan a la lista utilizando el operador +=.
    # La salida será una lista individuo que contiene un elemento aleatorio de 
    # genetic_pool[0] y otros elementos seleccionados aleatoriamente de 
    # genetic_pool[1].
    
    individuo = []
    individuo += [np.random.choice(genetic_pool[0])]
    individuo += list(np.random.choice(genetic_pool[1],ind_size-1))
    #print(individuo)
    decimal = listToDecimal(individuo)
    
    poblacion = []

    for i in range(100):
        individuo = []
        individuo += [np.random.choice(genetic_pool[0])]
        individuo += list(np.random.choice(genetic_pool[1],ind_size-1))
        poblacion.append(individuo)
    poblacion[:10]
    #print(poblacion)
    
    # Finalmente, genero una población llena de individuos con genes aleatorios.
    # De esta población vamos a elegir los mejores para reproducirlos.
    # Abajo observarás cómo se encuentra repartida la población.

    # for individuo in poblacion:
    #     x = listToDecimal(individuo)
    #     y = fx(x)
    #     plt.plot(x,y,'x')
    # plt.plot(x_axis,y_axis)
    # plt.show()
    
    
    # Fitness
    
    # Aquí vamos a medir el éxito del individuo para cumplir con el objetivo 
    # y determinar la probabilidad que tendrá de reproducirse.
    # Ya que queremos maximizar una función, aquellos individuos que produzcan
    # un valor más alto en y serán seleccionados como los mejores.
    fitness =[]
    
    #extraigo los valores de y para medir su éxito
    for individuo in poblacion:
        x = listToDecimal(individuo)
        y = fx(x)
        fitness += [y]
    
    #convierto fitnees en un vector para realizar operaciones más fácilmente
    fitness = np.array(fitness)

    #divido todos los valores de y para la suma total y así obtener valores entre 0 y 1
    fitness=fitness/fitness.sum()
    
    #print(fitness)
    
    # Todos los valores de y se dividen entre la suma total para obtener una probabilidad. 
    # La operación realizada arriba sirve para representar porcentajes. Los números más 
    # grandes producen un porcentaje mayor,y ya que la probabilidad se mide entre 0 y 1, 
    # esto ya noos permite darle una probabilidad mayor de reproducirse a los mejores 
    # individuos, aquellos que obtuvieron números más grandes en y.

    # Entrecruzamiento
    
    # consiste en mezclar los genes de los mejores individuos. En este caso vamos a elegir 
    # dos padres al azar, de acuerdo a la probabilidad del fitness, para que produzcan un 
    # nuevo individuo.
    
    # Es importante no descartar de lleno a los peores individuos ya que, tal vez, su código
    # genético puede servir para lograr algo mejor más tarde.
    
    # Luego, se elige un cross_point a partir del cual se van a combinar los genes de los 
    # padres. Se va a copiar los genes del primer padre hasta este punto, y luego los genes 
    # del otro padre.
    
    size_poblacion = len(poblacion)
    #hijos
    offspring = []
    for i in range(size_poblacion//2):
        parents = np.random.choice(size_poblacion, 2, p=fitness)
        cross_point = np.random.randint(ind_size)
        offspring += [poblacion[parents[0]][:cross_point] + poblacion[parents[1]][cross_point:]]
        offspring += [poblacion[parents[1]][:cross_point] + poblacion[parents[0]][cross_point:]]


    # Generar hijos en un algoritmo genético. La variable poblacion contiene la población 
    # actual de individuos, y se desea crear una lista offspring que contendrá los hijos 
    # generados a partir de los padres seleccionados aleatoriamente. 
    
    # Se utiliza un bucle for para generar una cantidad de hijos igual a la mitad del tamaño
    # de la población actual (size_poblacion // 2). Dentro del bucle, se seleccionan aleatoriamente 
    # dos padres de la población utilizando np.random.choice con el parámetro p que indica la 
    # probabilidad de selección basada en la aptitud de cada individuo (fitness). Luego, se 
    # selecciona aleatoriamente un punto de cruce (cross_point) utilizando np.random.randint en el 
    # rango de índices válidos para los cromosomas. Finalmente, se crea un hijo combinando los 
    # cromosomas de los padres en el punto de cruce y se agrega a la lista offspring. 
    
    # La salida será una lista offspring que contiene los hijos generados mediante el cruce
    # de los padres seleccionados aleatoriamente.
    
    # print(offspring[:10])
    poblacion = offspring
    # for individuo in poblacion:
    #     x = listToDecimal(individuo)
    #     y = fx(x)
    #     plt.plot(x,y,'x')
    # plt.plot(x_axis,y_axis)
    # plt.show()
    
    # En una sola generación podemos ver que los individuos ya se están concentrando en puntos más 
    # altos de la función. A medida que pasen las generaciones el objetivo es que todos los individuos 
    # convergan alrededor de 1.
    
    # Mutaciones
    # Otro aspecto interesante del algoritmo genético es la probabilidad de que se produzca una mutación 
    # en los individuos. Lo que significa que el ADN de cualquier individuo puede producir un gen que 
    # no viene de sus padres.
    
    # El individuo [0,5,3,8,9] puede mutar a [0,5,4,8,9] y esto mantiene la puerta abierta hacia el mejoramiento 
    # de la población.
    
    mutate(poblacion,0.005,genetic_pool)
    poblacion[:10]
    
    for individuo in poblacion:
        x = listToDecimal(individuo)
        y = fx(x)
        plt.plot(x,y,'x')
    plt.plot(x_axis,y_axis)
    
    # No hay un gran cambio debido a que la probabilidad de mutar es baja, pero así es como tiene que mantenerse, 
    # de lo contrario sería imposible llegar a una solución. La mutación debe ser muy poco frecuente.
    
    # Resultado final
    generaciones = 100

    for _ in range(generaciones):
    
        fitness =[]

    #extraigo los valores de y para medir su éxito
        for individuo in poblacion:
            x = listToDecimal(individuo)
            y = fx(x)
            fitness += [y]

    #convierto fitnees en un vector para realizar operaciones
    #más fácilmente
    fitness = np.array(fitness)

    #divido todos los valores de y para la suma total
    #y así obtener valores entre 0 y 1
    fitness=fitness/fitness.sum()    
        
    
    # se reproducen los mejores individuos
    offspring = []
    for i in range(size_poblacion//2):
        parents = np.random.choice(size_poblacion, 2, p=fitness)
        cross_point = np.random.randint(ind_size)
        offspring += [poblacion[parents[0]][:cross_point] + poblacion[parents[1]][cross_point:]]
        offspring += [poblacion[parents[1]][:cross_point] + poblacion[parents[0]][cross_point:]]    
    
    poblacion = offspring
    
    #####
    # MUTACIONES
    ####
     
    mutate(poblacion,0.005,genetic_pool)
    
    for individuo in poblacion:
        x = listToDecimal(individuo)
        y = fx(x)
        plt.plot(x,y,'x')
    plt.plot(x_axis,y_axis)
    
    # Después de 100 generaciones podemos ver que la mayoría de los individuos se concentran 
    # en el valor máximo. 
    # Y si vemos el mejor individuo vamos a ver que su valor es cercano a 1.
    
    np.where(fitness == fitness.max())
    listToDecimal(poblacion[41])
    fx(listToDecimal(poblacion[41]))
    
    #Podríamos dejarlos reproducirse por unas cuántas generaciones más a ver si siguen avanzando hacia el máximo.
    #Después de 300 generaciones ya casi todos los individuos están cerca del máximo.
    listToDecimal(poblacion[np.where(fitness == fitness.max())[0][0]]) 
    #El mejor individuo sigue siendo sólo cercano a 1, pues su ADN no es muy estable y es imposible obtener el valor 
    #exacto de 1. Sin embargo, esto es posible con un ADN binario.
    
if __name__ == '__main__':
    main()