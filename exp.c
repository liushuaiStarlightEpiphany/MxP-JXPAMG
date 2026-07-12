#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dirent.h>
#include <time.h>

#define PI 3.141592653589793

static int sgn(double x) { return (x>0)-(x<0); }
static int cmp_d(const void*a,const void*b){double d=*(const double*)a-*(const double*)b;return(d>0)-(d<0);}

void algo2(const double m[6],double l[3],double eps){
    double m0=fmin(fmin(fabs(m[0]),fabs(m[1])),fabs(m[2]));
    double r1=2*fabs(m[3])/m0,r2=2*fabs(m[4])/m0,r3=2*fabs(m[5])/m0;
    if(r1+r2+r3<eps){l[0]=m[0];l[1]=m[1];l[2]=m[2];qsort(l,3,sizeof(double),cmp_d);return;}
    if(r1+r2<eps){l[0]=m[0];double d=sqrt((m[1]-m[2])*(m[1]-m[2])+4*m[5]*m[5]);l[1]=(m[1]+m[2]+d)/2;l[2]=(m[1]+m[2]-d)/2;qsort(l,3,sizeof(double),cmp_d);return;}
    if(r1+r3<eps){l[1]=m[1];double d=sqrt((m[0]-m[2])*(m[0]-m[2])+4*m[4]*m[4]);l[0]=(m[0]+m[2]+d)/2;l[2]=(m[0]+m[2]-d)/2;qsort(l,3,sizeof(double),cmp_d);return;}
    if(r2+r3<eps){l[2]=m[2];double d=sqrt((m[0]-m[1])*(m[0]-m[1])+4*m[3]*m[3]);l[0]=(m[0]+m[1]+d)/2;l[1]=(m[0]+m[1]-d)/2;qsort(l,3,sizeof(double),cmp_d);return;}
    double t3=(m[0]+m[1]+m[2])/3,p=m[0]*m[1]+m[0]*m[2]+m[1]*m[2]-m[3]*m[3]-m[4]*m[4]-m[5]*m[5]-3*t3*t3;
    double q=m[0]*m[1]*m[2]+2*m[3]*m[5]*m[4]-m[0]*m[5]*m[5]-m[1]*m[4]*m[4]-m[2]*m[3]*m[3]+2*t3*t3*t3-t3*(m[0]*m[1]+m[0]*m[2]+m[1]*m[2]-m[3]*m[3]-m[4]*m[4]-m[5]*m[5]);
    double dp=p/3,dq=-q/2,dis=dq*dq+dp*dp*dp;
    if(dis>=0){double s=sqrt(dis),u=cbrt(dq+s),v=cbrt(dq-s);l[0]=u+v+t3;l[1]=-(u+v)/2+t3;l[2]=l[1];l[1]+=sqrt(3)/2*(u-v);l[2]-=sqrt(3)/2*(u-v);}
    else{double phi=acos(dq/sqrt(-dp*dp*dp))/3,r=2*sqrt(-dp);l[0]=r*cos(phi)+t3;l[1]=r*cos(phi+2*PI/3)+t3;l[2]=r*cos(phi+4*PI/3)+t3;}
    qsort(l,3,sizeof(double),cmp_d);
}

int get_case(double m[6],double eps){
    double m0=fmin(fmin(fabs(m[0]),fabs(m[1])),fabs(m[2]));
    double r1=2*fabs(m[3])/m0,r2=2*fabs(m[4])/m0,r3=2*fabs(m[5])/m0;
    if(r1+r2+r3<eps)return 1;if(r1+r2<eps)return 2;if(r1+r3<eps)return 3;if(r2+r3<eps)return 4;return 5;
}

char outpath[1024];

void write_csv(const char *group, const char *epsname, int cas, double err) {
    FILE *fp = fopen(outpath, "a");
    if(fp){fprintf(fp,"%s,%s,%d,%.15e\n",group,epsname,cas,err);fclose(fp);}
}

int read_mtx(const char *fname, double m[6]){
    FILE *fp=fopen(fname,"r");if(!fp)return 0;
    for(int i=0;i<6;i++)if(fscanf(fp,"%lf",&m[i])!=1){fclose(fp);return 0;}
    fclose(fp);return 1;
}

void scan_dir(const char *path, double *mats, int *n, int max_n) {
    DIR *dp=opendir(path);if(!dp)return;
    struct dirent *entry;char fpath[1024];
    while((entry=readdir(dp))&&*n<max_n){
        if(!strstr(entry->d_name,".txt"))continue;
        sprintf(fpath,"%s/%s",path,entry->d_name);
        double m[6];
        if(read_mtx(fpath,m)){for(int i=0;i<6;i++)mats[(*n)*6+i]=m[i];(*n)++;}
    }
    closedir(dp);
}

int main(){
    const char *groups[]={"diag","ps1","ps2","ps3","general","dual",
                          "sweep/2.5e-3","sweep/2.5e-4","sweep/2.5e-5","sweep/2.5e-6","sweep/2.5e-7",NULL};
    char home[256];sprintf(home,"%s",getenv("HOME"));
    double epslist[]={1e-14,1e-12,1e-10,1e-8,1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2};
    char epsnames[][10]={"1e-14","1e-12","1e-10","1e-8","1e-7","5e-7","1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2"};

    // Init CSV
    sprintf(outpath,"%s/l_s/eigenvalue-3/result/exp_data.csv",home);
    FILE *fp=fopen(outpath,"w");
    if(fp){fprintf(fp,"group,eps,case,relerr\n");fclose(fp);}

    // Init result files
    sprintf(outpath,"%s/l_s/eigenvalue-3/result/case_distribution.txt",home);
    fp=fopen(outpath,"w");if(fp){fprintf(fp,"group,eps,C1,C2,C3,C4,C5,reducible,total,pct\n");fclose(fp);}
    sprintf(outpath,"%s/l_s/eigenvalue-3/result/error_stats.txt",home);
    fp=fopen(outpath,"w");if(fp){fprintf(fp,"group,eps,avg_err,max_err\n");fclose(fp);}
    sprintf(outpath,"%s/l_s/eigenvalue-3/result/timing.txt",home);
    fp=fopen(outpath,"w");if(fp){fprintf(fp,"group,algo2_ns,cardano_ns,ratio\n");fclose(fp);}

    for(int gi=0;groups[gi];gi++){
        char path[1024];
        sprintf(path,"%s/l_s/eigenvalue-3/matrix-test/%s",home,groups[gi]);
        int n=0;
        double *mats=(double*)malloc(10000*6*sizeof(double));
        scan_dir(path,mats,&n,10000);
        if(!n){free(mats);continue;}

        printf("\n========== %s (n=%d) ==========\n",groups[gi],n);

        for(int ei=0;ei<16;ei++){
            double eps=epslist[ei];
            int cases[6]={0};double sum_err=0,max_err=0;int valid=0;
            double l[3],lr[3];

            for(int i=0;i<n;i++){
                double *m=&mats[i*6];
                int cas=get_case(m,eps);cases[cas]++;
                algo2(m,l,eps);
                algo2(m,lr,0);
                double rel=0;
                for(int k=0;k<3;k++){
                    double e=fabs(l[k]-lr[k])/(fabs(lr[k])+1e-300);
                    if(e>rel)rel=e;
                }
                sum_err+=rel;if(rel>max_err)max_err=rel;valid++;

                sprintf(outpath,"%s/l_s/eigenvalue-3/result/exp_data.csv",home);
                write_csv(groups[gi],epsnames[ei],cas,rel);
            }

            int red=cases[1]+cases[2]+cases[3]+cases[4];
            printf("eps=%-5s: C1=%3d C2=%3d C3=%3d C4=%3d C5=%3d | red=%3d/%d(%2.0f%%) | avg_err=%.2e max_err=%.2e\n",
                   epsnames[ei],cases[1],cases[2],cases[3],cases[4],cases[5],
                   red,valid,100.0*red/valid,sum_err/valid,max_err);

            sprintf(outpath,"%s/l_s/eigenvalue-3/result/case_distribution.txt",home);
            fp=fopen(outpath,"a");
            if(fp){fprintf(fp,"%s,%s,%d,%d,%d,%d,%d,%d,%d,%.0f\n",
                        groups[gi],epsnames[ei],cases[1],cases[2],cases[3],cases[4],cases[5],
                        red,valid,100.0*red/valid);fclose(fp);}

            sprintf(outpath,"%s/l_s/eigenvalue-3/result/error_stats.txt",home);
            fp=fopen(outpath,"a");
            if(fp){fprintf(fp,"%s,%s,%.2e,%.2e\n",groups[gi],epsnames[ei],sum_err/valid,max_err);fclose(fp);}
        }

        // Timing
        int rep=20000;clock_t ca,cb;
        double ta=0,tc=0;double l[3];
        for(int i=0;i<n;i++){
            double *m=&mats[i*6];
            ca=clock();for(int r=0;r<rep;r++)algo2(m,l,1e-10);cb=clock();ta+=(double)(cb-ca);
            ca=clock();for(int r=0;r<rep;r++)algo2(m,l,0);cb=clock();tc+=(double)(cb-ca);
        }
        double ta_ns=ta/n/rep/CLOCKS_PER_SEC*1e9;
        double tc_ns=tc/n/rep/CLOCKS_PER_SEC*1e9;
        printf("  Timing: Algo2=%.1fns Cardano=%.1fns ratio=%.2f\n",ta_ns,tc_ns,tc_ns/ta_ns);

        sprintf(outpath,"%s/l_s/eigenvalue-3/result/timing.txt",home);
        fp=fopen(outpath,"a");
        if(fp){fprintf(fp,"%s,%.1f,%.1f,%.2f\n",groups[gi],ta_ns,tc_ns,tc_ns/ta_ns);fclose(fp);}

        free(mats);
    }

    printf("\nResults saved to ~/l_s/eigenvalue-3/result/\n");
    return 0;
}
