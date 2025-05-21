#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define OUTPUT_DIR "D:/INAOE/Doctorado/STM32/"  // Ajusta esta ruta según tu entorno

// Función que calcula las ecuaciones del sistema caótico de Lorenz
void lorenz_frac(const double x[3], double xdot[3]) {
    double sigma = 10.0;
    double beta = 8.0 / 3.0;
    double rho = 28.0;
    xdot[0] = sigma * (x[1] - x[0]);
    xdot[1] = rho * x[0] - x[1] - x[0] * x[2];
    xdot[2] = -beta * x[2] + x[0] * x[1];
}

// Función para el sistema HO2 (no se utiliza en este ejemplo)
void ho2_system(const double *x, double *xdot) {
    double a = 0.2, b = 0.1;
    // Asumiendo que x tiene 4 componentes
    xdot[0] = x[1] * x[2];
    xdot[1] = x[0] - x[1] - a * x[3];
    xdot[2] = 1 - x[0] * x[0];
    xdot[3] = b * x[1];
}

// Función para calcular los coeficientes binomiales con truncamiento decimal
// Se calcula un arreglo de longitud (mm + 1) según:
//   c[0] = 1
//   c[j] = c[j-1] * (1 - (1+alpha)/j),  para j = 1,..., mm
// Luego se truncan los valores usando factor = 10^decimal.
// Además, se escribe el arreglo en "lorenz_coeficientes.rnd" en OUTPUT_DIR.
double* binomial_coef(double alpha, int mm, int decimal, int *length) {
    int len = mm + 1;  // coeficientes de 0 a mm
    *length = len;
    double *c = (double*)malloc(len * sizeof(double));
    if (!c) {
        return NULL;
    }
    c[0] = 1.0;
    for (int j = 1; j < len; j++) {
        c[j] = c[j - 1] * (1.0 - (1.0 + alpha) / j);
    }
    double factor = pow(10, decimal);
    for (int j = 0; j < len; j++) {
        c[j] = trunc(c[j] * factor) / factor;
    }
    
    char file_name[256];
    snprintf(file_name, sizeof(file_name), "%slorenz_coeficientes.rnd", OUTPUT_DIR);
    FILE *fp = fopen(file_name, "w");
    if (fp) {
        for (int j = 0; j < len; j++) {
            fprintf(fp, "%.15f\n", c[j]);
        }
        fclose(fp);
    }
    return c;
}

// Función que implementa el método de Grünwald–Letnikov para la derivada fraccionaria.
// Se asume que la matriz de estados "x" está almacenada en forma de vector de tamaño (k+1)*d,
// donde cada fila tiene d elementos. Asimismo, x_t es el vector de tiempos (k+1 elementos).
// Parámetros:
//   h       : paso de integración
//   h_alpha : h^alpha
//   k       : número total de pasos
//   alpha   : orden fraccionario
//   nu      : parámetro (si nu != 1, se suma desde j = 1 hasta min(nu, i)-1; si nu == 1, se suma desde j = 1 hasta i-1)
//   d       : dimensión del sistema (3 para Lorenz)
//   mm      : parámetro para el cálculo de coeficientes (se calculan mm+1 coeficientes)
//   decimal : precisión decimal
//   m       : tamaño de la ventana de memoria (m = Lm/h)
// Para cada iteración i (1 <= i <= k) se realiza:
//   sum_x = sum_{j in idx} c[j] * x[i - j, :]
//   x[i, :] = lorenz_frac(x[i-1, :]) * h_alpha - sum_x
// Además, se almacenan los resultados en dos archivos.
double* grunwald_letnikov(double *x, double h, double h_alpha, int k, double alpha,
                           double *x_t, int nu, int d, int mm, int decimal, int m) {
    int coef_length = 0;
    double *c = binomial_coef(alpha, mm, decimal, &coef_length);
    if (!c) {
        return NULL;
    }
    
    char file_name1[256], file_name2[256];
    snprintf(file_name1, sizeof(file_name1), "%slorenz_variables.rnd", OUTPUT_DIR);
    snprintf(file_name2, sizeof(file_name2), "%slorenz_sumatoria.rnd", OUTPUT_DIR);
    
    FILE *fp1 = fopen(file_name1, "w");
    FILE *fp2 = fopen(file_name2, "w");
    if (!fp1 || !fp2) {
        free(c);
        if (fp1) fclose(fp1);
        if (fp2) fclose(fp2);
        return NULL;
    }
    
    // Bucle principal: i de 1 a k (se asume que x[0] ya contiene la condición inicial)
    for (int i = 1; i <= k; i++) {
        double sum_x[10]; // Se asume que d <= 10; en nuestro caso d == 3.
        for (int di = 0; di < d; di++) {
            sum_x[di] = 0.0;
        }
        
        // Determinar el rango de índices para la suma.
        // En el código Python:
        //   if i > 1:
        //       idx = range(1, min(nu, i))   si nu != 1, o range(1, i) si nu == 1
        int start = 1;
        int end;
        if (i > 1) {
            if (nu != 1) {
                end = (i < nu ? i : nu);  // min(nu, i)
            } else {
                end = i;
            }
            for (int j = start; j < end; j++) {
                for (int di = 0; di < d; di++) {
                    // x[(i - j), di] se almacena en x[(i - j)*d + di]
                    sum_x[di] += c[j] * x[(i - j) * d + di];
                }
            }
        }
        
        // Calcular la derivada evaluando las ecuaciones de Lorenz en el estado anterior.
        double prev_state[3], deriv[3];
        for (int di = 0; di < d; di++) {
            prev_state[di] = x[(i - 1) * d + di];
        }
        lorenz_frac(prev_state, deriv);
        
        // Actualizar el estado: x[i, :] = lorenz_frac(x[i-1, :]) * h_alpha - sum_x
        for (int di = 0; di < d; di++) {
            x[i * d + di] = deriv[di] * h_alpha - sum_x[di];
        }
        
        // Guardar en el archivo los datos del paso i-1
        fprintf(fp1, "%.3f\t", x_t[i - 1]);
        for (int di = 0; di < d; di++) {
            fprintf(fp1, "%.10f\t", x[(i - 1) * d + di]);
        }
        fprintf(fp1, "\n");
        
        // Guardar la suma calculada (sum_x) en otro archivo
        for (int di = 0; di < d; di++) {
            fprintf(fp2, "%.15f\t", sum_x[di]);
        }
        fprintf(fp2, "\n");
        
        // Imprimir por pantalla cada 100 iteraciones (opcional)
        if (i % 100 == 0 && d == 3) {
            printf("%.3f\t%.10f\t%.10f\t%.10f\n", x_t[i], x[i * d + 0], x[i * d + 1], x[i * d + 2]);
        }
    }
    
    fclose(fp1);
    fclose(fp2);
    free(c);
    return x;
}

int main(void) {
    // Parámetros de integración y del método
    int decimal = 10;            // Precisión decimal deseada
    double alpha = 0.98;         // Orden fraccionario
    int d = 3;                 // Dimensión del sistema (3 para Lorenz)
    double x0[3] = {0.1, 0.1, 0.1};  // Condición inicial
    
    double t0 = 0.0, tf = 100.0; // Intervalo de tiempo
    double h = 0.01;           // Paso de integración
    double h_alpha = pow(h, alpha);
    
    // Longitud de memoria Lm y cálculo de m = Lm/h
    double Lm = 1.0;           // Aquí se usa Lm = 1.0 (según el código Python)
    int m = (int)(Lm / h);
    
    // Número total de pasos: k = (tf - t0)/h
    int k = (int)((tf - t0) / h);
    
    // Se determinan nu y mm según el código Python:
    //   if k < m, entonces nu = 1 y mm = k; de lo contrario, nu = m y mm = m.
    int nu, mm;
    if (k < m) {
        nu = 1;
        mm = k;
    } else {
        nu = m;
        mm = m;
    }
    
    // Reservar memoria para la matriz de estados x (k+1 filas, d columnas)
    double *x = (double*)calloc((k + 1) * d, sizeof(double));
    if (!x) {
        return -1;
    }
    
    // Reservar y calcular el vector de tiempos t (k+1 elementos)
    double *t_arr = (double*)malloc((k + 1) * sizeof(double));
    if (!t_arr) {
        free(x);
        return -1;
    }
    for (int i = 0; i <= k; i++) {
        t_arr[i] = t0 + i * h;
    }
    
    // Establecer la condición inicial en la primera fila de x
    for (int i = 0; i < d; i++) {
        x[i] = x0[i];
    }
    
    // Ejecutar el método de Grünwald–Letnikov (la función retorna x, aunque se modifica in situ)
    grunwald_letnikov(x, h, h_alpha, k, alpha, t_arr, nu, d, mm, decimal, m);
    
    // En este ejemplo no se implementa la graficación en C.
    
    free(x);
    free(t_arr);
    
    return 0;
}
