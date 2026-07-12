#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793

static int sc(double x) { return (x > 0) - (x < 0); }
static int cd(const void *a, const void *b) { double d = *(const double *)a - *(const double *)b; return (d > 0) - (d < 0); }

void robust(const double A[6], double lambda[3]);

static void qrf(const double m[6], double l[3]) {
    double A[9] = {m[0], m[3], m[4], m[3], m[1], m[5], m[4], m[5], m[2]};
    for (int it = 0; it < 50; it++) {
        double sh = A[8];
        for (int i = 0; i < 3; i++) A[i*3+i] -= sh;
        for (int i = 0; i < 2; i++) {
            double a = A[i*3+i], b = A[(i+1)*3+i];
            double r = hypot(a, b), c = a / r, s = -b / r;
            for (int j = i; j < 3; j++) { double x = A[i*3+j], y = A[(i+1)*3+j]; A[i*3+j] = c*x - s*y; A[(i+1)*3+j] = s*x + c*y; }
            for (int j = 0; j < 3; j++) { double x = A[j*3+i], y = A[j*3+i+1]; A[j*3+i] = c*x - s*y; A[j*3+i+1] = s*x + c*y; }
        }
        for (int i = 0; i < 3; i++) A[i*3+i] += sh;
    }
    l[0] = A[0]; l[1] = A[4]; l[2] = A[8];
    qsort(l, 3, sizeof(double), cd);
}

static void conv(const double *m, double *r) {
    r[0] = m[0]; r[1] = m[3]; r[2] = m[4];
    r[3] = m[1]; r[4] = m[5]; r[5] = m[2];
}

int main() {
    // Test matrices with different diagonal spans
    double tests[][6] = {
        // Small span [1, 100]
        {50, 80, 20, 1, 0.5, 0.3},
        // Medium span [1, 1e6]
        {1e6, 1e3, 1, 100, 10, 1},
        // Large span [1, 1e12] - similar to paper
        {1e12, 1e6, 1, 1e5, 1e3, 1},
    };
    int nt = 3;
    
    for (int ti = 0; ti < nt; ti++) {
        double *m = tests[ti];
        double span_max = fmax(fmax(m[0], m[1]), m[2]);
        double span_min = fmin(fmin(m[0], m[1]), m[2]);
        
        double l_qr[3], l_rob[3], rr[6];
        conv(m, rr);
        robust(rr, l_rob);
        qsort(l_rob, 3, sizeof(double), cd);
        
        qrf(m, l_qr);
        qsort(l_qr, 3, sizeof(double), cd);
        
        double max_rel = 0;
        for (int k = 0; k < 3; k++) {
            double e = fabs(l_qr[k] - l_rob[k]) / (fabs(l_rob[k]) + 1e-300);
            if (e > max_rel) max_rel = e;
        }
        
        printf("Span [%.0e, %.0e] (ratio=%.0e): QR max rel err = %.2e\n", 
               span_min, span_max, span_max/span_min, max_rel);
        
        printf("  Robust: %10.6e %10.6e %10.6e\n", l_rob[0], l_rob[1], l_rob[2]);
        printf("  QR:     %10.6e %10.6e %10.6e\n", l_qr[0], l_qr[1], l_qr[2]);
        printf("\n");
    }
    return 0;
}
