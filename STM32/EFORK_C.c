#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Parámetros del problema y constantes
#define N1 10000            // número de pasos
#define T 100.0             // tiempo total de simulación
#define ALPHA 0.9935        // orden fraccionario
#define LM 10.0             // memoria corta

// Parámetros del sistema de Lorenz
#define SIGMA 10.0
#define RHO   28.0
#define BETA  (8.0/3.0)


// Función que calcula el sistema de Lorenz dado (x, y, z)
void lorenz_system(double x, double y, double z, double *dx, double *dy, double *dz) {
    *dx = SIGMA * (y - x);
    *dy = x * (RHO - z) - y;
    *dz = x * y - BETA * z;
}

// Función para calcular el término de memoria (memoria fraccionaria)
// Se suma desde j = max(0, k - nu) hasta j = k-1.
// Nota: Para eficiencia, se puede precomputar gamma(2 - ALPHA) y, si ALPHA es fijo, usar
// un “ring buffer” para guardar solo los últimos 'nu+1' puntos.
double memory_fractional(int k, double t, const double *v_arr, const double *vtn, double h, double alpha, int nu) {
    double sum = 0.0;
    int start_idx = (k - nu > 0) ? (k - nu) : 0;
    double gamma_term = tgamma(2.0 - alpha); // Si es posible, precomputarlo fuera de la función.
    for (int j = start_idx; j < k; j++) {
        double t0 = vtn[j];
        double t1 = vtn[j+1];
        double v1 = pow(t - t0, 1.0 - alpha);
        double v2 = pow(t - t1, 1.0 - alpha);
        double diff = (v_arr[j+1] - v_arr[j]);
        sum += diff * (v1 - v2);
    }
    return sum / (h * gamma_term);
}

int main(void) {
    // Paso de integración y tamaño de memoria (en número de pasos)
    double h = T / N1;
    double ha = pow(h, ALPHA);
    int nu = (int)(LM / h);

    // Precalcular constantes (usando tgamma de math.h)
    double gamma1 = tgamma(1.0 + ALPHA);
    double gamma2 = tgamma(1.0 + 2 * ALPHA);
    double gamma3 = tgamma(1.0 + 3 * ALPHA);
    double c2 = pow(1.0 / (2.0 * tgamma(1.0 + ALPHA)), 1.0 / ALPHA);
    double c3 = pow(1.0 / (4.0 * tgamma(1.0 + ALPHA)), 1.0 / ALPHA);
    double a21 = 1.0 / (2.0 * gamma1 * gamma1);
    double a31 = (gamma1 * gamma1 * tgamma(2 * ALPHA + 1) + 2.0 * (tgamma(2 * ALPHA + 1) * tgamma(2 * ALPHA + 1)) - tgamma(3 * ALPHA + 1))
                 / (4.0 * gamma1 * gamma1 * (2.0 * (tgamma(2 * ALPHA + 1) * tgamma(2 * ALPHA + 1)) - tgamma(3 * ALPHA + 1)));
    double a32 = - tgamma(2 * ALPHA + 1) / (4.0 * (2.0 * (tgamma(2 * ALPHA + 1) * tgamma(2 * ALPHA + 1)) - tgamma(3 * ALPHA + 1)));
    
    double w1 = (8.0 * pow(gamma1, 3) * pow(tgamma(1.0 + 2 * ALPHA), 2) - 6.0 * pow(gamma1, 3) * tgamma(1.0 + 3 * ALPHA) + tgamma(1.0 + 2 * ALPHA) * tgamma(1.0 + 3 * ALPHA))
                / (gamma1 * tgamma(1.0 + 2 * ALPHA) * tgamma(1.0 + 3 * ALPHA));
    double w2 = 2.0 * (gamma1 * gamma1) * (4.0 * pow(tgamma(1.0 + 2 * ALPHA), 2) - tgamma(1.0 + 3 * ALPHA))
                / (tgamma(1.0 + 2 * ALPHA) * tgamma(1.0 + 3 * ALPHA));
    double w3 = -8.0 * (gamma1 * gamma1) * (2.0 * pow(tgamma(1.0 + 2 * ALPHA), 2) - tgamma(1.0 + 3 * ALPHA))
                / (tgamma(1.0 + 2 * ALPHA) * tgamma(1.0 + 3 * ALPHA));

    // Reservar memoria para almacenar la solución
    double *vtn = (double *)malloc((N1 + 1) * sizeof(double));
    double *vxn = (double *)malloc((N1 + 1) * sizeof(double));
    double *vyn = (double *)malloc((N1 + 1) * sizeof(double));
    double *vzn = (double *)malloc((N1 + 1) * sizeof(double));
    if (!vtn || !vxn || !vyn || !vzn) {
        // Manejo de error de asignación
        return -1;
    }

    // Condiciones iniciales
    vtn[0] = 0.0;
    vxn[0] = -10.0;
    vyn[0] = -10.0;
    vzn[0] = 30.0;

    // Primer paso (n = 0) sin incluir la memoria
    double x_n = vxn[0], y_n = vyn[0], z_n = vzn[0];
    double dx, dy, dz;
    lorenz_system(x_n, y_n, z_n, &dx, &dy, &dz);
    double K1x = ha * dx;
    double K1y = ha * dy;
    double K1z = ha * dz;
    
    // K2 evaluado en (x + a21*K1, y + a21*K1, z + a21*K1)
    lorenz_system(x_n + a21 * K1x, y_n + a21 * K1y, z_n + a21 * K1z, &dx, &dy, &dz);
    double K2x = ha * dx;
    double K2y = ha * dy;
    double K2z = ha * dz;
    
    // K3 evaluado en (x + a31*K2 + a32*K1, …)
    lorenz_system(x_n + a31 * K2x + a32 * K1x,
                  y_n + a31 * K2y + a32 * K1y,
                  z_n + a31 * K2z + a32 * K1z,
                  &dx, &dy, &dz);
    double K3x = ha * dx;
    double K3y = ha * dy;
    double K3z = ha * dz;

    double x_n1 = x_n + w1 * K1x + w2 * K2x + w3 * K3x;
    double y_n1 = y_n + w1 * K1y + w2 * K2y + w3 * K3y;
    double z_n1 = z_n + w1 * K1z + w2 * K2z + w3 * K3z;

    vtn[1] = h;
    vxn[1] = x_n1;
    vyn[1] = y_n1;
    vzn[1] = z_n1;

    // Actualizar el estado
    x_n = x_n1;
    y_n = y_n1;
    z_n = z_n1;

    // Bucle principal para n >= 1
    for (int n = 1; n < N1; n++) {
        double tn = n * h;

        // Calcular los términos de memoria para cada variable
        double mem_x = memory_fractional(n, tn, vxn, vtn, h, ALPHA, nu);
        double mem_y = memory_fractional(n, tn, vyn, vtn, h, ALPHA, nu);
        double mem_z = memory_fractional(n, tn, vzn, vtn, h, ALPHA, nu);

        // Evaluar el sistema de Lorenz con el término de memoria
        lorenz_system(x_n, y_n, z_n, &dx, &dy, &dz);
        double f1 = dx - mem_x;
        double f2 = dy - mem_y;
        double f3 = dz - mem_z;
        double K1x = ha * f1;
        double K1y = ha * f2;
        double K1z = ha * f3;

        // Para K2 y K3, se deben repetir evaluaciones similares en tiempos tn + c2*h y tn + c3*h,
        // usando los estados intermedios. Aquí se muestra una versión simplificada.
        lorenz_system(x_n + a21 * K1x, y_n + a21 * K1y, z_n + a21 * K1z, &dx, &dy, &dz);
        double K2x = ha * dx;
        double K2y = ha * dy;
        double K2z = ha * dz;

        lorenz_system(x_n + a31 * K2x + a32 * K1x,
                      y_n + a31 * K2y + a32 * K1y,
                      z_n + a31 * K2z + a32 * K1z,
                      &dx, &dy, &dz);
        double K3x = ha * dx;
        double K3y = ha * dy;
        double K3z = ha * dz;

        // Actualizar el estado con la combinación de K1, K2 y K3
        x_n1 = x_n + w1 * K1x + w2 * K2x + w3 * K3x;
        y_n1 = y_n + w1 * K1y + w2 * K2y + w3 * K3y;
        z_n1 = z_n + w1 * K1z + w2 * K2z + w3 * K3z;
        
        double tn1 = (n + 1) * h;
        vtn[n + 1] = tn1;
        vxn[n + 1] = x_n1;
        vyn[n + 1] = y_n1;
        vzn[n + 1] = z_n1;

        // Actualizar para el siguiente paso
        x_n = x_n1;
        y_n = y_n1;
        z_n = z_n1;
    }

    // Guardar los resultados en un archivo de texto.
    // En un STM32 deberás reemplazar esta parte por las rutinas de FatFS o la librería de E/S correspondiente.

    FILE *fp = fopen("D:/INAOE/Doctorado/STM32/EFORK_c.txt", "w");

    if (fp != NULL) {

        for (int i = 0; i <= N1; i++) {
            fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", vtn[i], vxn[i], vyn[i], vzn[i]);
        }
        fclose(fp);
    } else {
        // Manejo de error en la apertura del archivo.
    }

    // Liberar la memoria asignada
    free(vtn);
    free(vxn);
    free(vyn);
    free(vzn);

    return 0;
}
