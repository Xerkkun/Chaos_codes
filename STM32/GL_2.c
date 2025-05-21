#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

char OUTPUT_DIR[256];

/* Tamaños máximos (ajustar según necesidades del sistema). */
#define MAX_M 10000   // Tamaño máximo de la ventana de memoria
#define D 3           // Dimensión del sistema (3 para Lorenz)

/*-----------------------------------------------------------------------------
 * Función que calcula las ecuaciones del sistema caótico de Lorenz
 * Entrada : x[3]  (estado actual)
 * Salida  : xdot[3] (derivada en ese punto)
 *---------------------------------------------------------------------------*/
void lorenz_frac(const double x[3], double xdot[3]) {
    double sigma = 10.0;
    double beta  = 8.0 / 3.0;
    double rho   = 28.0;
    xdot[0] = sigma * (x[1] - x[0]);
    xdot[1] = rho * x[0] - x[1] - x[0] * x[2];
    xdot[2] = -beta * x[2] + x[0] * x[1];
}

/*-----------------------------------------------------------------------------
 * Función para calcular los coeficientes binomiales (sin truncamiento)
 * Se calculan para índices j = 0..mm.
 * c[0] = 1
 * c[j] = c[j-1]*(1 - (1+alpha)/j), j=1..mm
 *---------------------------------------------------------------------------*/
void binomial_coef(double alpha, int mm, double *c) {
    c[0] = 1.0;
    for (int j = 1; j <= mm; j++) {
        c[j] = c[j - 1] * (1.0 - (1.0 + alpha) / j);
    }
}

/*-----------------------------------------------------------------------------
 * Método de Grünwald–Letnikov con principio de memoria corta
 * usando un buffer circular de tamaño (m+1) para almacenar estados.
 *
 * Parámetros:
 *   alpha, h, h_alpha  -> orden fraccionario, paso, paso^alpha
 *   m                  -> longitud de memoria (Lm/h)
 *   k                  -> número total de pasos
 *   c[]                -> coeficientes binomiales [0..m]
 *   x0[]               -> estado inicial (tamaño D)
 *   t0, tf, h          -> tiempo inicial, final, paso
 *---------------------------------------------------------------------------*/
void grunwald_letnikov(double alpha, double h, double h_alpha,
                       int m, int k, const double *c,
                       const double x0[], double t0, double tf, int d)
{
    /* Archivo de salida para las variables. */
    char fname[256];
    snprintf(fname, sizeof(fname), "%slorenz_variables.rnd", OUTPUT_DIR);
    FILE *fp = fopen(fname, "w");
    if (!fp) {
        fprintf(stderr, "No se pudo abrir el archivo de salida.\n");
        return;
    }

    /* Buffer circular para almacenar los últimos m+1 estados. */
    static double xbuf[MAX_M+1][D];
    /* Vector para almacenar la derivada en cada paso. */
    double deriv[D];
    /* Vector para la suma interna. */
    double sum_x[D];

    /* Inicializamos el buffer circular con el estado inicial en la posición 0. */
    for (int di = 0; di < d; di++) {
        xbuf[0][di] = x0[di];
    }

    /* Escribimos el estado inicial en el archivo. */
    double t = t0;
    fprintf(fp, "%.3f\t%.10f\t%.10f\t%.10f\n",
            t, xbuf[0][0], xbuf[0][1], xbuf[0][2]);

    /* Iteramos desde i=1 hasta k. i es el índice de paso. */
    for (int i = 1; i <= k; i++) {
        t = t0 + i*h;

        /* Índice circular para almacenar el nuevo estado. */
        int current_idx = i % (m+1);

        /* Inicializamos sum_x en 0. */
        for (int di = 0; di < d; di++) {
            sum_x[di] = 0.0;
        }

        /* Determinamos cuántos pasos atrás sumamos.
           Si i < m, sumamos hasta i-1; si i >= m, sumamos m términos. */
        int pasos_atras = (i < m) ? i : m;

        /* Realizamos la suma con los coeficientes c[j], j=1..pasos_atras. */
        for (int j = 1; j <= pasos_atras; j++) {
            /* Queremos acceder al estado i-j, pero en el buffer circular.
               i-j mod (m+1) -> índice en el buffer. */
            int idx_old = ( (i - j) ) % (m+1);

            for (int di = 0; di < d; di++) {
                sum_x[di] += c[j] * xbuf[idx_old][di];
            }
        }

        /* Calculamos la derivada en el estado anterior (i-1).
           El estado anterior en el buffer es (i-1) mod (m+1). */
        int prev_idx = (i - 1) % (m+1);
        lorenz_frac(xbuf[prev_idx], deriv);

        /* Actualizamos el nuevo estado con la fórmula:
             x[i] = deriv*h_alpha - sum_x
           En el buffer circular, esto se almacena en current_idx. */
        for (int di = 0; di < d; di++) {
            xbuf[current_idx][di] = deriv[di]*h_alpha - sum_x[di];
        }

        /* Escribimos en el archivo el estado anterior (i-1). */
        fprintf(fp, "%.3f\t%.10f\t%.10f\t%.10f\n",
                t, xbuf[prev_idx][0], xbuf[prev_idx][1], xbuf[prev_idx][2]);

        /* (Opcional) Imprimir cada 100 iteraciones. */
        if (i % 100 == 0) {
            printf("%.3f\t%.10f\t%.10f\t%.10f\n",
                   t, xbuf[current_idx][0], xbuf[current_idx][1], xbuf[current_idx][2]);
        }
    }

    fclose(fp);
}

/*-----------------------------------------------------------------------------
 * main
 *---------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
    /* Configurar OUTPUT_DIR desde la línea de comandos, si se proporciona. */
    if (argc > 1) {
        strncpy(OUTPUT_DIR, argv[1], sizeof(OUTPUT_DIR));
        OUTPUT_DIR[sizeof(OUTPUT_DIR)-1] = '\0';
    } else {
        strcpy(OUTPUT_DIR, "D:/INAOE/Doctorado/STM32/");
    }

    /* Parámetros de la simulación. */
    double alpha = 0.98;  // Orden fraccionario
    double t0 = 0.0;      // Tiempo inicial
    double tf = 500.0;    // Tiempo final
    double h = 0.01;      // Paso
    double h_alpha = pow(h, alpha);

    /* Longitud de memoria. */
    double Lm = 0.1;
    int m = (int)(Lm / h);
    if (m > MAX_M) {
        fprintf(stderr, "Error: m=%d excede MAX_M=%d.\n", m, MAX_M);
        return -1;
    }

    /* Número total de pasos. */
    int k = (int)((tf - t0)/h);

    /* Estado inicial. */
    double x0[D] = {0.1, 0.1, 0.1};

    int mm = m;

    static double c[MAX_M+1];
    binomial_coef(alpha, mm, c);

    grunwald_letnikov(alpha, h, h_alpha, m, k, c, x0, t0, tf, D);

    return 0;
}
