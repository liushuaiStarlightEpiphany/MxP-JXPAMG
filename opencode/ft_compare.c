#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dirent.h>
#include <time.h>

#define PI 3.141592653589793

static int sc(double x) { return (x > 0) - (x < 0); }

static int cd(const void *a, const void *b)
{
    double d = *(const double *)a - *(const double *)b;
    return (d > 0) - (d < 0);
}

void robust(const double A[6], double lambda[3]);

/* ==================== Algo2 ==================== */
static void a2(const double m[6], double l[3], double eps)
{
    double m0 = fmin(fmin(fabs(m[0]), fabs(m[1])), fabs(m[2]));
    double r1 = 2 * fabs(m[3]) / m0;
    double r2 = 2 * fabs(m[4]) / m0;
    double r3 = 2 * fabs(m[5]) / m0;

    if (r1 + r2 + r3 < eps) {
        l[0] = m[0]; l[1] = m[1]; l[2] = m[2];
        qsort(l, 3, sizeof(double), cd); return;
    }
    if (r1 + r2 < eps) {
        l[0] = m[0];
        double d = sqrt((m[1] - m[2]) * (m[1] - m[2]) + 4 * m[5] * m[5]);
        l[1] = (m[1] + m[2] + d) / 2; l[2] = (m[1] + m[2] - d) / 2;
        qsort(l, 3, sizeof(double), cd); return;
    }
    if (r1 + r3 < eps) {
        l[1] = m[1];
        double d = sqrt((m[0] - m[2]) * (m[0] - m[2]) + 4 * m[4] * m[4]);
        l[0] = (m[0] + m[2] + d) / 2; l[2] = (m[0] + m[2] - d) / 2;
        qsort(l, 3, sizeof(double), cd); return;
    }
    if (r2 + r3 < eps) {
        l[2] = m[2];
        double d = sqrt((m[0] - m[1]) * (m[0] - m[1]) + 4 * m[3] * m[3]);
        l[0] = (m[0] + m[1] + d) / 2; l[1] = (m[0] + m[1] - d) / 2;
        qsort(l, 3, sizeof(double), cd); return;
    }

    double rr[6];
    rr[0] = m[0]; rr[1] = m[3]; rr[2] = m[4];
    rr[3] = m[1]; rr[4] = m[5]; rr[5] = m[2];
    robust(rr, l);
    qsort(l, 3, sizeof(double), cd);
}

static void cardano_a2(const double m[6], double l[3])
{
    double t3 = (m[0] + m[1] + m[2]) / 3;
    double p = m[0]*m[1] + m[0]*m[2] + m[1]*m[2]
             - m[3]*m[3] - m[4]*m[4] - m[5]*m[5] - 3*t3*t3;
    double q = m[0]*m[1]*m[2] + 2*m[3]*m[5]*m[4]
             - m[0]*m[5]*m[5] - m[1]*m[4]*m[4] - m[2]*m[3]*m[3]
             + 2*t3*t3*t3 - t3*(m[0]*m[1] + m[0]*m[2] + m[1]*m[2]
             - m[3]*m[3] - m[4]*m[4] - m[5]*m[5]);
    double dp = p / 3, dq = -q / 2, dis = dq*dq + dp*dp*dp;
    if (dis >= 0) {
        double s = sqrt(dis), u = cbrt(dq + s), v = cbrt(dq - s);
        l[0] = u + v + t3; l[1] = -(u + v) / 2 + t3; l[2] = l[1];
        l[1] += sqrt(3) / 2 * (u - v); l[2] -= sqrt(3) / 2 * (u - v);
    } else {
        double phi = acos(dq / sqrt(-dp*dp*dp)) / 3, r = 2 * sqrt(-dp);
        l[0] = r * cos(phi) + t3;
        l[1] = r * cos(phi + 2*PI/3) + t3;
        l[2] = r * cos(phi + 4*PI/3) + t3;
    }
    qsort(l, 3, sizeof(double), cd);
}
static void jacobi(const double m[6], double l[3])
{
    double A[9] = {m[0], m[3], m[4], m[3], m[1], m[5], m[4], m[5], m[2]};
    for (int it = 0; it < 50; it++) {
        double mx = 0; int p = 0, q = 1;
        for (int i = 0; i < 3; i++)
            for (int j = i + 1; j < 3; j++)
                if (fabs(A[i*3+j]) > mx) { mx = fabs(A[i*3+j]); p = i; q = j; }
        if (mx < 1e-30) break;
        double th = (A[q*3+q] - A[p*3+p]) / (2 * A[p*3+q]);
        double t = sc(th) / (fabs(th) + sqrt(1 + th*th));
        double c = 1 / sqrt(1 + t*t), s = c * t, tau = s / (1 + c);
        double ap = A[p*3+p], aq = A[q*3+q], apq = A[p*3+q];
        A[p*3+p] = ap - t * apq; A[q*3+q] = aq + t * apq;
        A[p*3+q] = A[q*3+p] = 0;
        for (int r = 0; r < 3; r++)
            if (r != p && r != q) {
                double arp = A[r*3+p], arq = A[r*3+q];
                A[r*3+p] = A[p*3+r] = arp - s * (arq + tau * arp);
                A[r*3+q] = A[q*3+r] = arq + s * (arp - tau * arq);
            }
    }
    l[0] = A[0]; l[1] = A[4]; l[2] = A[8];
    qsort(l, 3, sizeof(double), cd);
}

/* ==================== Symmetric QL (tred2 + tqli from Numerical Recipes) ==================== */
#define SQR(a) ((a)*(a))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double pythag(double a, double b) {
    double absa = fabs(a), absb = fabs(b);
    if (absa > absb) return absa * sqrt(1.0 + SQR(absb / absa));
    else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

static void tred2(double *d, double *e, double *a) {
    int n = 3;
    double scale, hh, h, g, f;
    for (int i = n - 1; i > 0; i--) {
        int l = i - 1;
        h = scale = 0.0;
        if (l > 0) {
            for (int k = 0; k <= l; k++) scale += fabs(a[i*3+k]);
            if (scale == 0.0) e[i] = a[i*3+l];
            else {
                for (int k = 0; k <= l; k++) { a[i*3+k] /= scale; h += a[i*3+k] * a[i*3+k]; }
                f = a[i*3+l]; g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i] = scale * g; h -= f * g; a[i*3+l] = f - g; f = 0.0;
                for (int j = 0; j <= l; j++) {
                    a[j*3+i] = a[i*3+j] / h; g = 0.0;
                    for (int k = 0; k <= j; k++) g += a[j*3+k] * a[i*3+k];
                    for (int k = j + 1; k <= l; k++) g += a[k*3+j] * a[i*3+k];
                    e[j] = g / h; f += e[j] * a[i*3+j];
                }
                hh = f / (h + h);
                for (int j = 0; j <= l; j++) {
                    f = a[i*3+j]; e[j] = g = e[j] - hh * f;
                    for (int k = 0; k <= j; k++) a[j*3+k] -= (f * e[k] + g * a[i*3+k]);
                }
            }
        } else e[i] = a[i*3+l];
        d[i] = h;
    }
    d[0] = 0.0; e[0] = 0.0;
    for (int i = 0; i < n; i++) {
        int l = i;
        if (d[i] != 0.0)
            for (int j = 0; j < l; j++) {
                g = 0.0; for (int k = 0; k < l; k++) g += a[i*3+k] * a[k*3+j];
                for (int k = 0; k < l; k++) a[k*3+j] -= g * a[k*3+i];
            }
        d[i] = a[i*3+i]; a[i*3+i] = 1.0;
        for (int j = 0; j < l; j++) a[j*3+i] = a[i*3+j] = 0.0;
    }
}

static void tqli(double *d, double *e, int n) {
    double s, r, p, g, f, dd, c, b;
    for (int i = 1; i < n; i++) e[i-1] = e[i];
    e[n-1] = 0.0;
    for (int l = 0; l < n; l++) {
        int iter = 0, m, i;
        do {
            for (m = l; m < n - 1; m++) { dd = fabs(d[m]) + fabs(d[m + 1]); if (fabs(e[m]) + dd == dd) break; }
            if (m != l) {
                if (iter++ == 100) break;
                g = (d[l + 1] - d[l]) / (2.0 * e[l]); r = pythag(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                s = c = 1.0; p = 0.0;
                for (i = m - 1; i >= l; i--) {
                    f = s * e[i]; b = c * e[i];
                    e[i + 1] = (r = pythag(f, g));
                    if (r == 0.0) { d[i + 1] -= p; e[m] = 0.0; break; }
                    s = f / r; c = g / r;
                    g = d[i + 1] - p; r = (d[i] - g) * s + 2.0 * c * b;
                    d[i + 1] = g + (p = s * r); g = c * r - b;
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p; e[l] = g; e[m] = 0.0;
            }
        } while (m != l);
    }
}

static void sym_ql(const double m[6], double l[3]) {
    double a[9] = {m[0], m[3], m[4], m[3], m[1], m[5], m[4], m[5], m[2]};
    double d[3], e[3];
    tred2(d, e, a);
    tqli(d, e, 3);
    l[0] = d[0]; l[1] = d[1]; l[2] = d[2];
    qsort(l, 3, sizeof(double), cd);
}

/* ==================== Helpers ==================== */
static int gc(const double m[6], double eps)
{
    double m0 = fmin(fmin(fabs(m[0]), fabs(m[1])), fabs(m[2]));
    double r1 = 2 * fabs(m[3]) / m0;
    double r2 = 2 * fabs(m[4]) / m0;
    double r3 = 2 * fabs(m[5]) / m0;
    if (r1 + r2 + r3 < eps) return 1;
    if (r1 + r2 < eps) return 2;
    if (r1 + r3 < eps) return 3;
    if (r2 + r3 < eps) return 4;
    return 5;
}

static void cv(const double *m, double *r)
{
    r[0] = m[0]; r[1] = m[3]; r[2] = m[4];
    r[3] = m[1]; r[4] = m[5]; r[5] = m[2];
}

static double relerr(double a, double b)
{
    return fabs(a - b) / (fabs(b) + 1e-300);
}

typedef struct { const char *g; const char *e[8]; int n; } TG;

int main()
{
    char home[256]; sprintf(home, "%s", getenv("HOME"));
    char base[1024]; sprintf(base, "%s/l_s/eigenvalue-3/result", home);
    char p[1024];

    sprintf(p, "%s/compare_timing.txt", base);
    FILE *f = fopen(p, "w");
    if (f) { fprintf(f, "group,algo2_ns,cardano_ns,jacobi_ns,qr_ns\n"); fclose(f); }

    sprintf(p, "%s/compare_accuracy.txt", base);
    f = fopen(p, "w");
    if (f) { fprintf(f, "group,algo2_err,cardano_err,jacobi_err,qr_err\n"); fclose(f); }

    TG ts[] = {
        {"diag",            {"1e-4"}, 1},
        {"ps1",             {"1e-4"}, 1},
        {"ps2",             {"1e-4"}, 1},
        {"ps3",             {"1e-4"}, 1},
        {"general",         {"1e-4"}, 1},
        {"dual",            {"1e-4","1e-6"}, 2},
        {"sweep/2.5e-3",    {"1e-2"}, 1},
        {"sweep/2.5e-4",    {"1e-3"}, 1},
        {"sweep/2.5e-5",    {"1e-4"}, 1},
        {"sweep/2.5e-6",    {"1e-5"}, 1},
        {"sweep/2.5e-7",    {"1e-6"}, 1},
        {NULL, {NULL}, 0}
    };

    for (int gi = 0; ts[gi].g; gi++) {
        sprintf(p, "%s/l_s/eigenvalue-3/matrix-test/%s", home, ts[gi].g);
        DIR *dp = opendir(p);
        if (!dp) continue;
        struct dirent *en; char fp2[1024];
        double *ms = malloc(10000 * 6 * sizeof(double));
        int n = 0;

        while ((en = readdir(dp)) && n < 10000) {
            if (!strstr(en->d_name, ".txt")) continue;
            sprintf(fp2, "%s/%s", p, en->d_name);
            FILE *f2 = fopen(fp2, "r");
            if (!f2) continue;
            double m[6]; int ok = 1;
            for (int i = 0; i < 6; i++)
                if (fscanf(f2, "%lf", &m[i]) != 1) ok = 0;
            fclose(f2);
            if (ok) {
                for (int i = 0; i < 6; i++) ms[n*6+i] = m[i];
                n++;
            }
        }
        closedir(dp);
        if (!n) { free(ms); continue; }

        printf("\n========== %s (n=%d) ==========\n", ts[gi].g, n);

        double sa2 = 0, sca = 0, sjac = 0, sqr = 0;
        double ea2 = 0, eca = 0, ejac = 0, eqr = 0;

        for (int ei = 0; ei < ts[gi].n; ei++) {
            double eps = atof(ts[gi].e[ei]);

            for (int i = 0; i < n; i++) {
                double *m = &ms[i*6];
                double l[3], lr[3], rr[6];
                cv(m, rr); robust(rr, lr);
                qsort(lr, 3, sizeof(double), cd);

                a2(m, l, eps);
                qsort(l, 3, sizeof(double), cd);
                double e = fmax(fmax(relerr(l[0],lr[0]), relerr(l[1],lr[1])), relerr(l[2],lr[2]));
                ea2 += e;

                cardano_a2(m, l); /* Cardano */
                qsort(l, 3, sizeof(double), cd);
                e = fmax(fmax(relerr(l[0],lr[0]), relerr(l[1],lr[1])), relerr(l[2],lr[2]));
                eca += e;

                jacobi(m, l);
                qsort(l, 3, sizeof(double), cd);
                e = fmax(fmax(relerr(l[0],lr[0]), relerr(l[1],lr[1])), relerr(l[2],lr[2]));
                ejac += e;

                sym_ql(m, l);
                qsort(l, 3, sizeof(double), cd);
                e = fmax(fmax(relerr(l[0],lr[0]), relerr(l[1],lr[1])), relerr(l[2],lr[2]));
                eqr += e;
            }

            double avga2 = ea2 / n, avgca = eca / n, avgjac = ejac / n, avgqr = eqr / n;
            printf("eps=%-5s: avg_err  a2=%.2e cardano=%.2e jacobi=%.2e qr=%.2e\n",
                   ts[gi].e[ei], avga2, avgca, avgjac, avgqr);

            sprintf(fp2, "%s/compare_accuracy.txt", base);
            f = fopen(fp2, "a");
            if (f) {
                fprintf(f, "%s_%s,%.2e,%.2e,%.2e,%.2e\n",
                        ts[gi].g, ts[gi].e[ei], avga2, avgca, avgjac, avgqr);
                fclose(f);
            }
        }

        /* Timing */
        int rep = 20000;
        clock_t ca, cb;
        double ta = 0, tc = 0, tj = 0, tq = 0;
        double l[3];

        for (int i = 0; i < n; i++) {
            double *m = &ms[i*6];

            ca = clock(); for (int r = 0; r < rep; r++) a2(m, l, 1e-10); cb = clock();
            ta += (double)(cb - ca);

            ca = clock(); for (int r = 0; r < rep; r++) cardano_a2(m, l); cb = clock();
            tc += (double)(cb - ca);

            ca = clock(); for (int r = 0; r < rep/4; r++) jacobi(m, l); cb = clock();
            tj += (double)(cb - ca);

            ca = clock(); for (int r = 0; r < rep/4; r++) sym_ql(m, l); cb = clock();
            tq += (double)(cb - ca);
        }

        double t_a2 = ta / n / rep / CLOCKS_PER_SEC * 1e9;
        double t_ca = tc / n / rep / CLOCKS_PER_SEC * 1e9;
        double t_ja = tj / n / (rep/4) / CLOCKS_PER_SEC * 1e9;
        double t_qr = tq / n / (rep/4) / CLOCKS_PER_SEC * 1e9;
        printf("  Timing: a2=%.1fns cardano=%.1fns jacobi=%.1fns qr=%.1fns\n",
               t_a2, t_ca, t_ja, t_qr);

        sprintf(fp2, "%s/compare_timing.txt", base);
        f = fopen(fp2, "a");
        if (f) {
            fprintf(f, "%s,%.1f,%.1f,%.1f,%.1f\n", ts[gi].g, t_a2, t_ca, t_ja, t_qr);
            fclose(f);
        }

        free(ms);
    }

    printf("\nDone. Results in result/compare_*.txt\n");
    return 0;
}
