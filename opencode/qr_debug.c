#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int sc(double x) { return (x > 0) - (x < 0); }
static int cd(const void *a, const void *b) { double d = *(const double *)a - *(const double *)b; return (d > 0) - (d < 0); }
void robust(const double A[6], double lambda[3]);

static void qrf_debug(const double m[6], double l[3], int print) {
    double A[9] = {m[0], m[3], m[4], m[3], m[1], m[5], m[4], m[5], m[2]};
    for (int it = 0; it < 50; it++) {
        double sh = A[8];
        for (int i = 0; i < 3; i++) A[i*3+i] -= sh;
        for (int i = 0; i < 2; i++) {
            double a = A[i*3+i], b = A[(i+1)*3+i];
            double r = hypot(a, b);
            if (r == 0) continue;
            double c = a / r, s = -b / r;
            for (int j = i; j < 3; j++) { double x = A[i*3+j], y = A[(i+1)*3+j]; A[i*3+j] = c*x - s*y; A[(i+1)*3+j] = s*x + c*y; }
            for (int j = 0; j < 3; j++) { double x = A[j*3+i], y = A[j*3+i+1]; A[j*3+i] = c*x - s*y; A[j*3+i+1] = s*x + c*y; }
        }
        for (int i = 0; i < 3; i++) A[i*3+i] += sh;
        
        if (print && it < 5) {
            printf("  it %d: diag = %.6e %.6e %.6e, off = %.2e %.2e %.2e\n", 
                   it, A[0], A[4], A[8], fabs(A[1]), fabs(A[2]), fabs(A[5]));
        }
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
    // Single test matrix: large span
    double m[6] = {1e12, 1e6, 1.0, 1e5, 100, 0.01};
    double rr[6], lr[3], lq[3];
    conv(m, rr);
    robust(rr, lr);
    qsort(lr, 3, sizeof(double), cd);
    
    printf("=== Single matrix test ===\n");
    printf("Matrix: a11=%.0e a22=%.0e a33=%.0e a12=%.0e a13=%.0e a23=%.0e\n", m[0], m[1], m[2], m[3], m[4], m[5]);
    printf("\nQR iterations:\n");
    qrf_debug(m, lq, 1);
    qsort(lq, 3, sizeof(double), cd);
    
    printf("\nRobust: %.15e  %.15e  %.15e\n", lr[0], lr[1], lr[2]);
    printf("QR:     %.15e  %.15e  %.15e\n", lq[0], lq[1], lq[2]);
    
    double e0 = relerr(lq[0], lr[0]);
    double e1 = relerr(lq[1], lr[1]);
    double e2 = relerr(lq[2], lr[2]);
    printf("Rel err per eigenvalue: %.2e  %.2e  %.2e\n", e0, e1, e2);
    printf("Max rel err: %.2e\n", fmax(fmax(e0, e1), e2));
    
    // Test with more iterations
    printf("\n=== With 200 QR iterations ===\n");
    {
        double A[9] = {m[0], m[3], m[4], m[3], m[1], m[5], m[4], m[5], m[2]};
        for (int it = 0; it < 200; it++) {
            double sh = A[8];
            for (int i = 0; i < 3; i++) A[i*3+i] -= sh;
            for (int i = 0; i < 2; i++) {
                double a = A[i*3+i], b = A[(i+1)*3+i];
                double r = hypot(a, b);
                if (r == 0) continue;
                double c = a / r, s = -b / r;
                for (int j = i; j < 3; j++) { double x = A[i*3+j], y = A[(i+1)*3+j]; A[i*3+j] = c*x - s*y; A[(i+1)*3+j] = s*x + c*y; }
                for (int j = 0; j < 3; j++) { double x = A[j*3+i], y = A[j*3+i+1]; A[j*3+i] = c*x - s*y; A[j*3+i+1] = s*x + c*y; }
            }
            for (int i = 0; i < 3; i++) A[i*3+i] += sh;
        }
        double l200[3] = {A[0], A[4], A[8]};
        qsort(l200, 3, sizeof(double), cd);
        printf("QR(200): %.15e  %.15e  %.15e\n", l200[0], l200[1], l200[2]);
        printf("Max rel err: %.2e\n", fmax(fmax(relerr(l200[0],lr[0]), relerr(l200[1],lr[1])), relerr(l200[2],lr[2])));
    }

    return 0;
}
