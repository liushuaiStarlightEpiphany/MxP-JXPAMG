#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dirent.h>
#include <time.h>

static int sc(double x)
{
    return (x > 0) - (x < 0);
}

static int cd(const void *a, const void *b)
{
    double d = *(const double *)a - *(const double *)b;
    return (d > 0) - (d < 0);
}

void robust(const double A[6], double lambda[3]);

static void a2(const double m[6], double l[3], double eps)
{
    double m0 = fmin(fmin(fabs(m[0]), fabs(m[1])), fabs(m[2]));
    double r1 = 2 * fabs(m[3]) / m0;
    double r2 = 2 * fabs(m[4]) / m0;
    double r3 = 2 * fabs(m[5]) / m0;

    if (r1 + r2 + r3 < eps) {
        l[0] = m[0]; l[1] = m[1]; l[2] = m[2];
        qsort(l, 3, sizeof(double), cd);
        return;
    }
    if (r1 + r2 < eps) {
        l[0] = m[0];
        double d = sqrt((m[1] - m[2]) * (m[1] - m[2]) + 4 * m[5] * m[5]);
        l[1] = (m[1] + m[2] + d) / 2;
        l[2] = (m[1] + m[2] - d) / 2;
        qsort(l, 3, sizeof(double), cd);
        return;
    }
    if (r1 + r3 < eps) {
        l[1] = m[1];
        double d = sqrt((m[0] - m[2]) * (m[0] - m[2]) + 4 * m[4] * m[4]);
        l[0] = (m[0] + m[2] + d) / 2;
        l[2] = (m[0] + m[2] - d) / 2;
        qsort(l, 3, sizeof(double), cd);
        return;
    }
    if (r2 + r3 < eps) {
        l[2] = m[2];
        double d = sqrt((m[0] - m[1]) * (m[0] - m[1]) + 4 * m[3] * m[3]);
        l[0] = (m[0] + m[1] + d) / 2;
        l[1] = (m[0] + m[1] - d) / 2;
        qsort(l, 3, sizeof(double), cd);
        return;
    }

    double rr[6];
    rr[0] = m[0]; rr[1] = m[3]; rr[2] = m[4];
    rr[3] = m[1]; rr[4] = m[5]; rr[5] = m[2];
    robust(rr, l);
    qsort(l, 3, sizeof(double), cd);
}

static void cardano_a2(const double m[6], double l[3])
{
    double t3  = (m[0] + m[1] + m[2]) / 3;
    double p   = m[0]*m[1] + m[0]*m[2] + m[1]*m[2]
               - m[3]*m[3] - m[4]*m[4] - m[5]*m[5] - 3*t3*t3;
    double q   = m[0]*m[1]*m[2] + 2*m[3]*m[5]*m[4]
               - m[0]*m[5]*m[5] - m[1]*m[4]*m[4] - m[2]*m[3]*m[3]
               + 2*t3*t3*t3 - t3*(m[0]*m[1] + m[0]*m[2] + m[1]*m[2]
               - m[3]*m[3] - m[4]*m[4] - m[5]*m[5]);
    double dp  = p / 3;
    double dq  = -q / 2;
    double dis = dq*dq + dp*dp*dp;

    if (dis >= 0) {
        double s = sqrt(dis);
        double u = cbrt(dq + s);
        double v = cbrt(dq - s);
        l[0] = u + v + t3;
        l[1] = -(u + v) / 2 + t3;
        l[2] = l[1];
        l[1] += sqrt(3) / 2 * (u - v);
        l[2] -= sqrt(3) / 2 * (u - v);
    } else {
        double phi = acos(dq / sqrt(-dp*dp*dp)) / 3;
        double r   = 2 * sqrt(-dp);
        l[0] = r * cos(phi) + t3;
        l[1] = r * cos(phi + 2*M_PI/3) + t3;
        l[2] = r * cos(phi + 4*M_PI/3) + t3;
    }
    qsort(l, 3, sizeof(double), cd);
}

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

typedef struct {
    const char *g;
    const char *e[8];
    int n;
} TG;

int main()
{
    char home[256];
    sprintf(home, "%s", getenv("HOME"));
    char base[1024];
    sprintf(base, "%s/l_s/eigenvalue-3/result", home);
    char p[1024];

    sprintf(p, "%s/exp_data.csv", base);
    FILE *f = fopen(p, "w");
    if (f) { fprintf(f, "group,eps,case,relerr\n"); fclose(f); }

    sprintf(p, "%s/case_distribution.txt", base);
    f = fopen(p, "w");
    if (f) { fprintf(f, "group,eps,C1,C2,C3,C4,C5,reducible,total,pct\n"); fclose(f); }

    sprintf(p, "%s/error_stats.txt", base);
    f = fopen(p, "w");
    if (f) { fprintf(f, "group,eps,avg_err,max_err\n"); fclose(f); }

    sprintf(p, "%s/timing.txt", base);
    f = fopen(p, "w");
    if (f) { fprintf(f, "group,algo2_ns,cardano_ns\n"); fclose(f); }

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

        struct dirent *en;
        char fp2[1024];
        double *ms = malloc(10000 * 6 * sizeof(double));
        int n = 0;

        while ((en = readdir(dp)) && n < 10000) {
            if (!strstr(en->d_name, ".txt")) continue;
            sprintf(fp2, "%s/%s", p, en->d_name);
            FILE *f2 = fopen(fp2, "r");
            if (!f2) continue;
            double m[6];
            int ok = 1;
            for (int i = 0; i < 6; i++)
                if (fscanf(f2, "%lf", &m[i]) != 1) ok = 0;
            fclose(f2);
            if (ok) {
                for (int i = 0; i < 6; i++) ms[n * 6 + i] = m[i];
                n++;
            }
        }
        closedir(dp);

        if (!n) { free(ms); continue; }

        printf("\n========== %s (n=%d) ==========\n", ts[gi].g, n);

        for (int ei = 0; ei < ts[gi].n; ei++) {
            double eps = atof(ts[gi].e[ei]);
            int cs[6] = {0};
            double se = 0, me = 0;
            int v = 0;
            double l[3], lr[3], rr[6];

            for (int i = 0; i < n; i++) {
                double *m = &ms[i * 6];
                int ca = gc(m, eps);
                cs[ca]++;
                a2(m, l, eps);
                cv(m, rr);
                robust(rr, lr);
                qsort(lr, 3, sizeof(double), cd);

                double rel = 0;
                for (int k = 0; k < 3; k++) {
                    double e = fabs(l[k] - lr[k]) / (fabs(lr[k]) + 1e-300);
                    if (e > rel) rel = e;
                }
                se += rel;
                if (rel > me) me = rel;
                v++;

                sprintf(fp2, "%s/exp_data.csv", base);
                f = fopen(fp2, "a");
                if (f) {
                    fprintf(f, "%s,%s,%d,%.15e\n", ts[gi].g, ts[gi].e[ei], ca, rel);
                    fclose(f);
                }
            }

            int red = cs[1] + cs[2] + cs[3] + cs[4];
            printf("eps=%-5s: C1=%3d C2=%3d C3=%3d C4=%3d C5=%3d | "
                   "red=%3d/%d(%2.0f%%) | avg_err=%.2e max_err=%.2e\n",
                   ts[gi].e[ei], cs[1], cs[2], cs[3], cs[4], cs[5],
                   red, v, 100.0 * red / v, se / v, me);

            sprintf(fp2, "%s/case_distribution.txt", base);
            f = fopen(fp2, "a");
            if (f) {
                fprintf(f, "%s,%s,%d,%d,%d,%d,%d,%d,%d,%.0f\n",
                        ts[gi].g, ts[gi].e[ei],
                        cs[1], cs[2], cs[3], cs[4], cs[5],
                        red, v, 100.0 * red / v);
                fclose(f);
            }

            sprintf(fp2, "%s/error_stats.txt", base);
            f = fopen(fp2, "a");
            if (f) {
                fprintf(f, "%s,%s,%.2e,%.2e\n", ts[gi].g, ts[gi].e[ei], se / v, me);
                fclose(f);
            }
        }

        int rep = 20000;
        clock_t ca, cb;
        double ta = 0, tc = 0;
        double l[3];

        for (int i = 0; i < n; i++) {
            double *m = &ms[i * 6];
            ca = clock();
            for (int r = 0; r < rep; r++) a2(m, l, 1e-10);
            cb = clock();
            ta += (double)(cb - ca);

            ca = clock();
            for (int r = 0; r < rep; r++) cardano_a2(m, l);
            cb = clock();
            tc += (double)(cb - ca);
        }

        double ta_ns = ta / n / rep / CLOCKS_PER_SEC * 1e9;
        double tc_ns = tc / n / rep / CLOCKS_PER_SEC * 1e9;
        printf("  Timing: Algo2=%.1fns Cardano=%.1fns\n", ta_ns, tc_ns);

        sprintf(fp2, "%s/timing.txt", base);
        f = fopen(fp2, "a");
        if (f) { fprintf(f, "%s,%.1f,%.1f\n", ts[gi].g, ta_ns, tc_ns); fclose(f); }

        free(ms);
    }

    printf("\nDone\n");
    return 0;
}
