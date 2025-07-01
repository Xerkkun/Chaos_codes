#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Parámetros de los sistemas
#define SIGMA 10.0
#define RHO   28.0
#define BETA  (8.0/3.0)

#define U 7.5
#define V 1.0
#define W 5.0

#define A 0.2
#define B 0.2
#define C 5.7

// Funciones para cada sistema
void lorenz_system(double x, double y, double z, double *dx, double *dy, double *dz) {
    *dx = SIGMA * (y - x);
    *dy = x * (RHO - z) - y;
    *dz = x * y - BETA * z;
}
void chen_system(double x, double y, double z, double *dx, double *dy, double *dz) {
    *dx = U * (y - x);
    *dy = (W - U) * x - x * z + W * y;
    *dz = x * y - V * z;
}
void rossler_system(double x, double y, double z, double *dx, double *dy, double *dz) {
    *dx = -y - z;
    *dy = x + A * y;
    *dz = B + z * (x - C);
}

// Término de memoria fraccional
double memory_fractional(int k, double t, const double *v_arr, const double *vtn, double h, double alpha, int nu) {
    double sum = 0.0;
    int start_idx = (k - nu > 0) ? (k - nu) : 0;
    double gamma_term = tgamma(2.0 - alpha);
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
    int opcion;
    char sistema[20] = "";
    double x0, y0, z0;
    double alpha, Lm, t_f, h;
    int N1, nu;

    // Menú de selección
    printf("Seleccione el sistema caotico a simular:\n");
    printf("1) Lorenz\n");
    printf("2) Rossler\n");
    printf("3) Chen\n");
    printf("Opcion (1/2/3): ");
    scanf("%d", &opcion);

    // Asignar función según selección
    void (*system_func)(double, double, double, double *, double *, double *);
    switch (opcion) {
        case 1: strcpy(sistema, "lorenz");   system_func = lorenz_system; break;
        case 2: strcpy(sistema, "rossler");  system_func = rossler_system; break;
        case 3: strcpy(sistema, "chen");     system_func = chen_system; break;
        default: printf("Seleccion invalida.\n"); return 1;
    }

    // Captura interactiva de parámetros
    printf("Ingrese orden fraccionario alpha: ");
    scanf("%lf", &alpha);

    printf("Ingrese paso h: ");
    scanf("%lf", &h);

    printf("Ingrese longitud de memoria Lm: ");
    scanf("%lf", &Lm);

    printf("Ingrese tiempo final t_f: ");
    scanf("%lf", &t_f);

    printf("Ingrese condiciones iniciales x0 y0 z0 (ejemplo: 0.1 0.1 0.1): ");
    scanf("%lf %lf %lf", &x0, &y0, &z0);

    // Calculo de parámetros derivados
    N1 = (int)ceil(t_f / h);
    nu = (int)(Lm / h);
    double ha = pow(h, alpha);

    // Constantes EFORK
    double gamma1 = tgamma(1.0 + alpha);
    double gamma2 = tgamma(1.0 + 2 * alpha);
    double gamma3 = tgamma(1.0 + 3 * alpha);
    double c2 = pow(1.0 / (2.0 * gamma1), 1.0 / alpha);
    double c3 = pow(1.0 / (4.0 * gamma1), 1.0 / alpha);
    double denom = 2.0 * gamma1 * gamma1 * gamma2 * gamma2 - gamma3 * gamma1 * gamma1;
    double a21 = 1.0 / (2.0 * gamma1 * gamma1);
    double a31 = (gamma1 * gamma1 * tgamma(2 * alpha + 1) + 2.0 * tgamma(2 * alpha + 1) * tgamma(2 * alpha + 1) - tgamma(3 * alpha + 1))
                 / (4.0 * gamma1 * gamma1 * (2.0 * tgamma(2 * alpha + 1) * tgamma(2 * alpha + 1) - tgamma(3 * alpha + 1)));
    double a32 = - tgamma(2 * alpha + 1) / (4.0 * (2.0 * tgamma(2 * alpha + 1) * tgamma(2 * alpha + 1) - tgamma(3 * alpha + 1)));
    double w1 = (8.0 * pow(gamma1, 3) * pow(tgamma(1.0 + 2 * alpha), 2) - 6.0 * pow(gamma1, 3) * tgamma(1.0 + 3 * alpha) + tgamma(1.0 + 2 * alpha) * tgamma(1.0 + 3 * alpha))
                / (gamma1 * tgamma(1.0 + 2 * alpha) * tgamma(1.0 + 3 * alpha));
    double w2 = 2.0 * gamma1 * gamma1 * (4.0 * pow(tgamma(1.0 + 2 * alpha), 2) - tgamma(1.0 + 3 * alpha))
                / (tgamma(1.0 + 2 * alpha) * tgamma(1.0 + 3 * alpha));
    double w3 = -8.0 * gamma1 * gamma1 * (2.0 * pow(tgamma(1.0 + 2 * alpha), 2) - tgamma(1.0 + 3 * alpha))
                / (tgamma(1.0 + 2 * alpha) * tgamma(1.0 + 3 * alpha));

    // Reservar memoria
    double *vtn = (double *)malloc((N1 + 1) * sizeof(double));
    double *vxn = (double *)malloc((N1 + 1) * sizeof(double));
    double *vyn = (double *)malloc((N1 + 1) * sizeof(double));
    double *vzn = (double *)malloc((N1 + 1) * sizeof(double));
    if (!vtn || !vxn || !vyn || !vzn) {
        printf("Error al asignar memoria.\n");
        return -1;
    }

    // Inicialización
    vtn[0] = 0.0;
    vxn[0] = x0;
    vyn[0] = y0;
    vzn[0] = z0;
    double x_n = x0, y_n = y0, z_n = z0;

    // Primer paso (n=0, sin memoria)
    double dx, dy, dz;
    system_func(x_n, y_n, z_n, &dx, &dy, &dz);
    double K1x = ha * dx, K1y = ha * dy, K1z = ha * dz;
    system_func(x_n + a21 * K1x, y_n + a21 * K1y, z_n + a21 * K1z, &dx, &dy, &dz);
    double K2x = ha * dx, K2y = ha * dy, K2z = ha * dz;
    system_func(x_n + a31 * K2x + a32 * K1x,
                y_n + a31 * K2y + a32 * K1y,
                z_n + a31 * K2z + a32 * K1z, &dx, &dy, &dz);
    double K3x = ha * dx, K3y = ha * dy, K3z = ha * dz;

    double x_n1 = x_n + w1 * K1x + w2 * K2x + w3 * K3x;
    double y_n1 = y_n + w1 * K1y + w2 * K2y + w3 * K3y;
    double z_n1 = z_n + w1 * K1z + w2 * K2z + w3 * K3z;
    vtn[1] = h;
    vxn[1] = x_n1;
    vyn[1] = y_n1;
    vzn[1] = z_n1;
    x_n = x_n1; y_n = y_n1; z_n = z_n1;

    // Bucle principal
    for (int n = 1; n < N1; n++) {
        double tn = n * h;

        // Término de memoria para cada variable
        double mem_x = memory_fractional(n, tn, vxn, vtn, h, alpha, nu);
        double mem_y = memory_fractional(n, tn, vyn, vtn, h, alpha, nu);
        double mem_z = memory_fractional(n, tn, vzn, vtn, h, alpha, nu);

        // Evaluar sistema con término de memoria
        system_func(x_n, y_n, z_n, &dx, &dy, &dz);
        double f1 = dx - mem_x, f2 = dy - mem_y, f3 = dz - mem_z;
        double K1x = ha * f1, K1y = ha * f2, K1z = ha * f3;

        system_func(x_n + a21 * K1x, y_n + a21 * K1y, z_n + a21 * K1z, &dx, &dy, &dz);
        double K2x = ha * dx, K2y = ha * dy, K2z = ha * dz;

        system_func(x_n + a31 * K2x + a32 * K1x,
                    y_n + a31 * K2y + a32 * K1y,
                    z_n + a31 * K2z + a32 * K1z, &dx, &dy, &dz);
        double K3x = ha * dx, K3y = ha * dy, K3z = ha * dz;

        // Actualizar estado
        x_n1 = x_n + w1 * K1x + w2 * K2x + w3 * K3x;
        y_n1 = y_n + w1 * K1y + w2 * K2y + w3 * K3y;
        z_n1 = z_n + w1 * K1z + w2 * K2z + w3 * K3z;

        double tn1 = (n + 1) * h;
        vtn[n + 1] = tn1;
        vxn[n + 1] = x_n1;
        vyn[n + 1] = y_n1;
        vzn[n + 1] = z_n1;

        x_n = x_n1; y_n = y_n1; z_n = z_n1;
    }

    // Guardar resultados en .txt
    char filename[128];
    snprintf(filename, sizeof(filename), "EFORK_%s_c.dat", sistema);
    FILE *fp = fopen(filename, "w");
    if (fp != NULL) {
        for (int i = 0; i <= N1; i++)
            fprintf(fp, "%.10f\t%.10f\t%.10f\t%.10f\n", vtn[i], vxn[i], vyn[i], vzn[i]);
        fclose(fp);
        printf("Archivo guardado: %s\n", filename);
    } else {
        printf("Error al guardar el archivo.\n");
    }


    free(vtn); free(vxn); free(vyn); free(vzn);
    return 0;
}
