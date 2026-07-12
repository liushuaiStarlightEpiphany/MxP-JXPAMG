#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>

#define PI 3.141592653589793

static int sgn(double x) { return (x > 0) - (x < 0); }
static int cmp_d(const void *a, const void *b) { double d = *(const double *)a - *(const double *)b; return (d > 0) - (d < 0); }

void algo2(const double m[6], double l[3], double eps) {
    double r1 = fabs(m[3]/m[0]) + fabs(m[3]/m[1]);
    double r2 = fabs(m[4]/m[0]) + fabs(m[4]/m[2]);
    double r3 = fabs(m[5]/m[1]) + fabs(m[5]/m[2]);
    if (r1+r2+r3 < eps) { l[0]=m[0]; l[1]=m[1]; l[2]=m[2]; qsort(l,3,sizeof(double),cmp_d); return; }
    if (r1+r2 < eps) { l[0]=m[0]; double d=sqrt((m[1]-m[2])*(m[1]-m[2])+4*m[5]*m[5]); l[1]=(m[1]+m[2]+d)/2; l[2]=(m[1]+m[2]-d)/2; qsort(l,3,sizeof(double),cmp_d); return; }
    if (r1+r3 < eps) { l[1]=m[1]; double d=sqrt((m[0]-m[2])*(m[0]-m[2])+4*m[4]*m[4]); l[0]=(m[0]+m[2]+d)/2; l[2]=(m[0]+m[2]-d)/2; qsort(l,3,sizeof(double),cmp_d); return; }
    if (r2+r3 < eps) { l[2]=m[2]; double d=sqrt((m[0]-m[1])*(m[0]-m[1])+4*m[3]*m[3]); l[0]=(m[0]+m[1]+d)/2; l[1]=(m[0]+m[1]-d)/2; qsort(l,3,sizeof(double),cmp_d); return; }
    double t3=(m[0]+m[1]+m[2])/3;
    double p=m[0]*m[1]+m[0]*m[2]+m[1]*m[2]-m[3]*m[3]-m[4]*m[4]-m[5]*m[5]-3*t3*t3;
    double q=m[0]*m[1]*m[2]+2*m[3]*m[5]*m[4]-m[0]*m[5]*m[5]-m[1]*m[4]*m[4]-m[2]*m[3]*m[3]+2*t3*t3*t3-t3*(m[0]*m[1]+m[0]*m[2]+m[1]*m[2]-m[3]*m[3]-m[4]*m[4]-m[5]*m[5]);
    double dp=p/3, dq=-q/2, dis=dq*dq+dp*dp*dp;
    if (dis>=0) { double s=sqrt(dis),u=cbrt(dq+s),v=cbrt(dq-s); l[0]=u+v+t3; l[1]=-(u+v)/2+t3; l[2]=l[1]; l[1]+=sqrt(3)/2*(u-v); l[2]-=sqrt(3)/2*(u-v); }
    else { double phi=acos(dq/sqrt(-dp*dp*dp))/3, r=2*sqrt(-dp); l[0]=r*cos(phi)+t3; l[1]=r*cos(phi+2*PI/3)+t3; l[2]=r*cos(phi+4*PI/3)+t3; }
    qsort(l,3,sizeof(double),cmp_d);
}

int get_case(double m[6], double eps) {
    double r1 = fabs(m[3]/m[0]) + fabs(m[3]/m[1]);
    double r2 = fabs(m[4]/m[0]) + fabs(m[4]/m[2]);
    double r3 = fabs(m[5]/m[1]) + fabs(m[5]/m[2]);
    if (r1+r2+r3 < eps) return 1;
    if (r1+r2 < eps) return 2;
    if (r1+r3 < eps) return 3;
    if (r2+r3 < eps) return 4;
    return 5;
}

int read_mtx(const char *fname, double m[6]) {
    FILE *fp = fopen(fname, "r");
    if (!fp) return 0;
    for (int i = 0; i < 6; i++)
        if (fscanf(fp, "%lf", &m[i]) != 1) { fclose(fp); return 0; }
    fclose(fp); return 1;
}

void scan_dir(const char *path, double *mats, int *n, int max_n) {
    DIR *dp = opendir(path);
    if (!dp) return;
    struct dirent *entry;
    char fpath[1024];
    while ((entry = readdir(dp)) && *n < max_n) {
        if (strstr(entry->d_name, ".txt") == NULL) continue;
        sprintf(fpath, "%s/%s", path, entry->d_name);
        double m[6];
        if (read_mtx(fpath, m)) {
            for (int i = 0; i < 6; i++) mats[(*n) * 6 + i] = m[i];
            (*n)++;
        }
    }
    closedir(dp);
}

int main() {
    const char *groups[] = {"diag", "ps1", "ps2", "ps3", "general", "dual", "sweep/1e-2", "sweep/1e-3", "sweep/1e-4", "sweep/1e-5", "sweep/1e-6", NULL};
    double eps_list[] = {1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2};
    
    for (int gi = 0; groups[gi]; gi++) {
        char base[1024];
        sprintf(base, "%s/%s", getenv("HOME"), "l_s/eigenvalue-3/matrix-test");
        char path[1024];
        sprintf(path, "%s/%s", base, groups[gi]);
        
        // Count files
        int n = 0;
        double *mats = malloc(10000 * 6 * sizeof(double));
        scan_dir(path, mats, &n, 10000);
        if (n == 0) { free(mats); continue; }
        
        printf("\n========== %s (n=%d) ==========\n", groups[gi], n);
        
        for (int ei = 0; ei < 7; ei++) {
            double eps = eps_list[ei];
            int cases[6] = {0};
            double sum_err = 0, max_err = 0;
            int valid = 0;
            double l[3];
            
            for (int i = 0; i < n; i++) {
                double *m = &mats[i * 6];
                int cas = get_case(m, eps);
                cases[cas]++;
                algo2(m, l, eps);
                // reference: algo2 with eps=0 (Cardano)
                double lr[3];
                algo2(m, lr, 0);
                double rel = 0;
                for (int k = 0; k < 3; k++) {
                    double e = fabs(l[k]-lr[k]) / (fabs(lr[k])+1e-300);
                    if (e > rel) rel = e;
                }
                sum_err += rel;
                if (rel > max_err) max_err = rel;
                valid++;
            }
            
            int red = cases[1]+cases[2]+cases[3]+cases[4];
            printf("eps=%.0e: C1=%d C2=%d C3=%d C4=%d C5=%d | reduce=%d/%d(%.0f%%) | avg_err=%.2e max_err=%.2e\n",
                   eps, cases[1], cases[2], cases[3], cases[4], cases[5],
                   red, valid, 100.0*red/valid,
                   sum_err/valid, max_err);
        }
        free(mats);
    }
    
    // Timing
    printf("\n========== TIMING ==========\n");
    for (int gi = 0; groups[gi]; gi++) {
        char path[1024];
        sprintf(path, "%s/%s/%s", getenv("HOME"), "l_s/eigenvalue-3/matrix-test", groups[gi]);
        int n = 0;
        double *mats = malloc(10000 * 6 * sizeof(double));
        scan_dir(path, mats, &n, 10000);
        if (n == 0) { free(mats); continue; }
        
        int rep = 20000;
        clock_t ca, cb;
        double ta = 0, tc = 0;
        double l[3];
        for (int i = 0; i < n; i++) {
            double *m = &mats[i*6];
            ca = clock(); for (int r = 0; r < rep; r++) algo2(m, l, 1e-10); cb = clock(); ta += (double)(cb-ca);
            ca = clock(); for (int r = 0; r < rep; r++) algo2(m, l, 0); cb = clock(); tc += (double)(cb-ca);
        }
        printf("%-12s Algo2=%.1f ns  Cardano=%.1f ns  ratio=%.2f\n",
               groups[gi],
               ta/n/rep/CLOCKS_PER_SEC*1e9,
               tc/n/rep/CLOCKS_PER_SEC*1e9,
               tc/ta);
        free(mats);
    }
    
    return 0;
}
