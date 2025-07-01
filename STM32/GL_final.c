#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_M 10000   // Tamaño máximo de la ventana de memoria
#define D 3           // Dimensión del sistema

/*-----------------------------------------------------------------------------
 * Funciones que calculan las ecuaciones de los sistemas caóticos (fraccionales)
 *---------------------------------------------------------------------------*/
void lorenz_frac(const double x[3], double xdot[3]) {
    double sigma = 10.0;
    double beta  = 8.0 / 3.0;
    double rho   = 28.0;
    xdot[0] = sigma * (x[1] - x[0]);
    xdot[1] = rho * x[0] - x[1] - x[0] * x[2];
    xdot[2] = -beta * x[2] + x[0] * x[1];
}

void chen_frac(const double x[3], double xdot[3]) {
    double u = 7.5;
    double v = 1.0;
    double w = 5.0;
    xdot[0] = u * (x[1] - x[0]);
    xdot[1] = (w - u) * x[0] - x[0] * x[2] + w * x[1];
    xdot[2] = x[0] * x[1] - v * x[2];
}

void rossler_frac(const double x[3], double xdot[3]) {
    double a = 0.2;
    double b = 0.2;
    double c = 5.7;
    xdot[0] = -x[1] - x[2];
    xdot[1] = x[0] + a * x[1];
    xdot[2] = b + x[2] * (x[0] - c);
}

/*-----------------------------------------------------------------------------
 * Función para calcular los coeficientes binomiales (sin truncamiento)
 *---------------------------------------------------------------------------*/
void binomial_coef(double alpha, int mm, double *c) {
    c[0] = 1.0;
    for (int j = 1; j <= mm; j++) {
        c[j] = c[j - 1] * (1.0 - (1.0 + alpha) / j);
    }
}

/*-----------------------------------------------------------------------------
 * Método de Grünwald–Letnikov con principio de memoria corta
 * y almacenamiento circular de estados.
 *---------------------------------------------------------------------------*/
void grunwald_letnikov(double alpha, double h, double h_alpha,
                       int m, int k, const double *c,
                       const double x0[], double t0, double tf, int d,
                       void (*system_frac)(const double[], double[]),
                       const char *sistema)
{
    // Nombre de archivo de salida según el sistema
    char filename[256];
    snprintf(filename, sizeof(filename), "GL_%s_c.dat", sistema);
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "No se pudo abrir el archivo de salida.\n");
        return;
    }

    static double xbuf[MAX_M+1][D];
    double deriv[D];
    double sum_x[D];

    // Inicializar buffer circular con estado inicial
    for (int di = 0; di < d; di++) {
        xbuf[0][di] = x0[di];
    }

    // Escribir estado inicial
    double t = t0;
    fprintf(fp, "%.10f\t%.10f\t%.10f\t%.10f\n",
            t, xbuf[0][0], xbuf[0][1], xbuf[0][2]);

    // Iterar pasos
    for (int i = 1; i <= k; i++) {
        t = t0 + i * h;
        int current_idx = i % (m + 1);

        // Sumatoria de memoria
        for (int di = 0; di < d; di++) sum_x[di] = 0.0;
        int pasos_atras = (i < m) ? i : m;
        for (int j = 1; j <= pasos_atras; j++) {
            int idx_old = (i - j) % (m + 1);
            for (int di = 0; di < d; di++) {
                sum_x[di] += c[j] * xbuf[idx_old][di];
            }
        }

        // Cálculo de derivadas
        int prev_idx = (i - 1) % (m + 1);
        system_frac(xbuf[prev_idx], deriv);

        // Actualizar estado
        for (int di = 0; di < d; di++) {
            xbuf[current_idx][di] = deriv[di] * h_alpha - sum_x[di];
        }

        // Escribir nuevo estado
        fprintf(fp, "%.10f\t%.10f\t%.10f\t%.10f\n",
                t, xbuf[current_idx][0], xbuf[current_idx][1], xbuf[current_idx][2]);
    }

    fclose(fp);
}

int main(void) {
    int opcion;
    char sistema[20] = "";
    double x0[D];
    double alpha, Lm, t_f, h;
    int m, k;

    void (*system_frac)(const double[], double[]);

    // Menú de selección de sistema
    printf("Seleccione el sistema caotico a simular:\n");
    printf("1) Lorenz\n");
    printf("2) Rossler\n");
    printf("3) Chen\n");
    printf("Opcion (1/2/3): ");
    if (scanf("%d", &opcion) != 1) {
        fprintf(stderr, "Entrada invalida.\n");
        return 1;
    }
    switch (opcion) {
        case 1: strcpy(sistema, "lorenz");   system_frac = lorenz_frac;  break;
        case 2: strcpy(sistema, "rossler");  system_frac = rossler_frac; break;
        case 3: strcpy(sistema, "chen");     system_frac = chen_frac;    break;
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
    scanf("%lf %lf %lf", &x0[0], &x0[1], &x0[2]);

    // Cálculo de parámetros derivados
    m = (int)(Lm / h);
    if (m > MAX_M) {
        fprintf(stderr, "Error: m=%d excede MAX_M=%d.\n", m, MAX_M);
        return -1;
    }
    k = (int)(t_f / h);
    double h_alpha = pow(h, alpha);

    static double c[MAX_M+1];
    binomial_coef(alpha, m, c);

    // Ejecutar GL
    grunwald_letnikov(alpha, h, h_alpha, m, k, c, x0, 0.0, t_f, D,
                       system_frac, sistema);

    return 0;
}
