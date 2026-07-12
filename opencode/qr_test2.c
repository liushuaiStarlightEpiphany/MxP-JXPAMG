#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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

double relerr(double a, double b) { return fabs(a - b) / (fabs(b) + 1e-300); }

int main() {
    srand48(time(NULL));
    int n = 500;
    double sum_qr = 0, max_qr = 0;
    double sum_ja = 0, max_ja = 0;
    int n_ok = 0;

    for (int i = 0; i < n; i++) {
        // Generate matrix similar to paper: diagonal in [1, 1e12] log-uniform
        double d[3];
        for (int k = 0; k < 3; k++)
            d[k] = pow(10, drand48() * 12);

        double m0 = fmin(fmin(d[0], d[1]), d[2]);
        // Off-diagonal entries: up to ~49% of min diagonal (semi-strong diagonal dominance)
        double a12 = (drand48() * 2 - 1) * 0.49 * m0;
        double a13 = (drand48() * 2 - 1) * 0.49 * m0;
        double a23 = (drand48() * 2 - 1) * 0.49 * m0;
        
        double m[6] = {d[0], d[1], d[2], a12, a13, a23};
        double rr[6], lr[3], lq[3];
        conv(m, rr);
        robust(rr, lr);
        qsort(lr, 3, sizeof(double), cd);

        qrf(m, lq);
        qsort(lq, 3, sizeof(double), cd);

        double e = fmax(fmax(relerr(lq[0],lr[0]), relerr(lq[1],lr[1])), relerr(lq[2],lr[2]));
        sum_qr += e;
        if (e > max_qr) max_qr = e;
        
        // Also compute condition number estimate: max/min eigenvalue
        double span_ratio = lr[2] / fmax(lr[0], 1e-300);
        if (span_ratio < 1e6) n_ok++;
    }

    printf("QR over %d random matrices (span 1~1e12):\n", n);
    printf("  Avg rel err: %.2e\n", sum_qr / n);
    printf("  Max rel err: %.2e\n", max_qr);
    
    // Test with smaller span for comparison
    n = 500;
    sum_qr = 0; max_qr = 0;
    for (int i = 0; i < n; i++) {
        double d[3];
        for (int k = 0; k < 3; k++)
            d[k] = pow(10, drand48() * 2); // span [1, 100]

        double m0 = fmin(fmin(d[0], d[1]), d[2]);
        double a12 = (drand48() * 2 - 1) * 0.49 * m0;
        double a13 = (drand48() * 2 - 1) * 0.49 * m0;
        double a23 = (drand48() * 2 - 1) * 0.49 * m0;
        
        double m[6] = {d[0], d[1], d[2], a12, a13, a23};
        double rr[6], lr[3], lq[3];
        conv(m, rr);
        robust(rr, lr);
        qsort(lr, 3, sizeof(double), cd);

        qrf(m, lq);
        qsort(lq, 3, sizeof(double), cd);

        double e = fmax(fmax(relerr(lq[0],lr[0]), relerr(lq[1],lr[1])), relerr(lq[2],lr[2]));
        sum_qr += e;
        if (e > max_qr) max_qr = e;
    }
    printf("\nQR over %d random matrices (span 1~100):\n", n);
    printf("  Avg rel err: %.2e\n", sum_qr / n);
    printf("  Max rel err: %.2e\n", max_qr);

    return 0;
}
