#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dirent.h>
#include <time.h>
#define PI acos(-1.0)
static int sc(double x){return(x>0)-(x<0);}
static int cd(const void*a,const void*b){double d=*(const double*)a-*(const double*)b;return(d>0)-(d<0);}

/* algo2 */
void a2(const double m[6],double l[3],double eps){
    double m0=fmin(fmin(fabs(m[0]),fabs(m[1])),fabs(m[2]));
    double r1=2*fabs(m[3])/m0,r2=2*fabs(m[4])/m0,r3=2*fabs(m[5])/m0;
    if(r1+r2+r3<eps){l[0]=m[0];l[1]=m[1];l[2]=m[2];qsort(l,3,sizeof(double),cd);return;}
    if(r1+r2<eps){l[0]=m[0];double d=sqrt((m[1]-m[2])*(m[1]-m[2])+4*m[5]*m[5]);l[1]=(m[1]+m[2]+d)/2;l[2]=(m[1]+m[2]-d)/2;qsort(l,3,sizeof(double),cd);return;}
    if(r1+r3<eps){l[1]=m[1];double d=sqrt((m[0]-m[2])*(m[0]-m[2])+4*m[4]*m[4]);l[0]=(m[0]+m[2]+d)/2;l[2]=(m[0]+m[2]-d)/2;qsort(l,3,sizeof(double),cd);return;}
    if(r2+r3<eps){l[2]=m[2];double d=sqrt((m[0]-m[1])*(m[0]-m[1])+4*m[3]*m[3]);l[0]=(m[0]+m[1]+d)/2;l[1]=(m[0]+m[1]-d)/2;qsort(l,3,sizeof(double),cd);return;}
    double t3=(m[0]+m[1]+m[2])/3,p=m[0]*m[1]+m[0]*m[2]+m[1]*m[2]-m[3]*m[3]-m[4]*m[4]-m[5]*m[5]-3*t3*t3;
    double q=m[0]*m[1]*m[2]+2*m[3]*m[5]*m[4]-m[0]*m[5]*m[5]-m[1]*m[4]*m[4]-m[2]*m[3]*m[3]+2*t3*t3*t3-t3*(m[0]*m[1]+m[0]*m[2]+m[1]*m[2]-m[3]*m[3]-m[4]*m[4]-m[5]*m[5]);
    double dp=p/3,dq=-q/2,dis=dq*dq+dp*dp*dp;
    if(dis>=0){double s=sqrt(dis),u=cbrt(dq+s),v=cbrt(dq-s);l[0]=u+v+t3;l[1]=-(u+v)/2+t3;l[2]=l[1];l[1]+=sqrt(3)/2*(u-v);l[2]-=sqrt(3)/2*(u-v);}
    else{double phi=acos(dq/sqrt(-dp*dp*dp))/3,r=2*sqrt(-dp);l[0]=r*cos(phi)+t3;l[1]=r*cos(phi+2*PI/3)+t3;l[2]=r*cos(phi+4*PI/3)+t3;}
    qsort(l,3,sizeof(double),cd);
}

/* robust helpers */
double tr_r(const double*A){return A[0]+A[3]+A[5];}
void sh_r(const double*A,double*B){for(int i=0;i<6;i++)B[i]=A[i];double a=tr_r(A)/3;B[0]-=a;B[3]-=a;B[5]-=a;}
double de_r(const double*A){return A[0]*A[3]*A[5]+A[1]*A[4]*A[2]+A[2]*A[1]*A[4]-A[2]*A[3]*A[2]-A[4]*A[4]*A[0]-A[1]*A[1]*A[5];}
void mv_r(const double*A,const double*v,double*r){r[0]=A[0]*v[0]+A[1]*v[1]+A[2]*v[2];r[1]=A[1]*v[0]+A[3]*v[1]+A[4]*v[2];r[2]=A[2]*v[0]+A[4]*v[1]+A[5]*v[2];}
void su_r(const double*A,const double*B,double*C){for(int i=0;i<6;i++)C[i]=A[i]-B[i];}
void sv_r(const double*a,const double*b,double*c){c[0]=a[0]-b[0];c[1]=a[1]-b[1];c[2]=a[2]-b[2];}
double mo_r(const double*v){return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);}
double do_r(const double*a,const double*b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
void sc_r(double n,double*v,double*r){r[0]=v[0]*n;r[1]=v[1]*n;r[2]=v[2]*n;}
int sg_r(double x){return(x>0)-(x<0);}

/* robust (Scherzinger) */
void robust(const double A[6],double*lambda){
    if(A[1]==0&&A[2]==0&&A[4]==0){lambda[0]=A[0];lambda[1]=A[3];lambda[2]=A[5];return;}
    double t3=tr_r(A)/3,A1[6];sh_r(A,A1);
    double J2=A1[1]*A1[1]+A1[2]*A1[2]+A1[4]*A1[4]-A1[0]*A1[3]-A1[0]*A1[5]-A1[3]*A1[5];
    double J3=de_r(A1);double alpha=10;
    if(J2!=0){double t=J3/(2.0*sqrt(J2*J2*J2/27.0));if(t>1)t=1;alpha=acos(t)/3.0;}
    double A1l1;
    if(alpha<PI/6)A1l1=2.0*sqrt(J2/3.0)*cos(alpha);
    else if(alpha>PI/6&&alpha!=10)A1l1=2.0*sqrt(J2/3.0)*cos(alpha+2*PI/3);
    if(alpha==PI/6||alpha==10){lambda[0]=2.0*sqrt(J2/3.0)*cos(alpha)+t3;lambda[1]=2.0*sqrt(J2/3.0)*cos(alpha+2*PI/3)+t3;lambda[2]=2.0*sqrt(J2/3.0)*cos(alpha+4*PI/3)+t3;return;}
    double e1[3]={1,0,0},e2[3]={0,1,0},e3[3]={0,0,1};
    double D[6]={A1l1,0,0,A1l1,0,A1l1},A3[6];su_r(A1,D,A3);
    double r1[3],r2[3],r3[3],vx[3];mv_r(A3,e1,r1);mv_r(A3,e2,r2);mv_r(A3,e3,r3);
    if(mo_r(r1)<mo_r(r2)){for(int i=0;i<3;i++){vx[i]=r1[i];r1[i]=r2[i];r2[i]=vx[i];}}
    if(mo_r(r1)<mo_r(r3)){for(int i=0;i<3;i++){vx[i]=r1[i];r1[i]=r3[i];r3[i]=vx[i];}}
    double s1[3];sc_r(1.0/mo_r(r1),r1,s1);
    double q=do_r(s1,r2);sc_r(q,s1,vx);double t2[3];sv_r(r2,vx,t2);
    q=do_r(s1,r3);sc_r(q,s1,vx);double t3v[3];sv_r(r3,vx,t3v);
    if(mo_r(t2)<mo_r(t3v)){for(int j=0;j<3;j++){vx[j]=t2[j];t2[j]=t3v[j];t3v[j]=vx[j];}}
    double s2[3];sc_r(1.0/mo_r(t2),t2,s2);
    double s11[3],s12[3];mv_r(A1,s1,s11);mv_r(A1,s2,s12);
    double A4[6]={A1l1,do_r(s1,s11),do_r(s1,s12),do_r(s2,s11),do_r(s2,s12)};
    double d=A4[1]-A4[4];
    double A1l2=(A4[1]+A4[4])/2.0-sg_r(A4[1]-A4[4])*sqrt(d*d+4*A4[2]*A4[3])/2.0;
    double A1l3=(A4[1]+A4[4])/2.0+sg_r(A4[1]-A4[4])*sqrt(d*d+4*A4[2]*A4[3])/2.0;
    lambda[0]=A1l1+t3;lambda[1]=A1l2+t3;lambda[2]=A1l3+t3;
}

int get_case(double m[6],double eps){
    double m0=fmin(fmin(fabs(m[0]),fabs(m[1])),fabs(m[2]));
    double r1=2*fabs(m[3])/m0,r2=2*fabs(m[4])/m0,r3=2*fabs(m[5])/m0;
    if(r1+r2+r3<eps)return 1;if(r1+r2<eps)return 2;if(r1+r3<eps)return 3;if(r2+r3<eps)return 4;return 5;
}

void conv(double*m,double*r){r[0]=m[0];r[1]=m[3];r[2]=m[4];r[3]=m[1];r[4]=m[5];r[5]=m[2];}

int main(){
    char home[256];sprintf(home,"%s",getenv("HOME"));
    char base[1024];sprintf(base,"%s/l_s/eigenvalue-3/result",home);
    char path[1024];
    sprintf(path,"%s/exp_data.csv",base);FILE*fp=fopen(path,"w");if(fp){fprintf(fp,"group,eps,case,relerr\n");fclose(fp);}
    sprintf(path,"%s/case_distribution.txt",base);fp=fopen(path,"w");if(fp){fprintf(fp,"group,eps,C1,C2,C3,C4,C5,reducible,total,pct\n");fclose(fp);}
    sprintf(path,"%s/error_stats.txt",base);fp=fopen(path,"w");if(fp){fprintf(fp,"group,eps,avg_err,max_err\n");fclose(fp);}
    sprintf(path,"%s/timing.txt",base);fp=fopen(path,"w");if(fp){fprintf(fp,"group,algo2_ns,cardano_ns,ratio\n");fclose(fp);}

    typedef struct{const char*g;const char*e[8];int n;}TG;
    TG ts[]={
        {"diag",{"1e-4"},1},{"ps1",{"1e-4"},1},{"ps2",{"1e-4"},1},{"ps3",{"1e-4"},1},{"general",{"1e-4"},1},
        {"dual",{"1e-4","1e-6"},2},{"sweep/2.5e-3",{"1e-2"},1},{"sweep/2.5e-4",{"1e-3"},1},{"sweep/2.5e-5",{"1e-4"},1},{"sweep/2.5e-6",{"1e-5"},1},{"sweep/2.5e-7",{"1e-6"},1},{NULL,{NULL},0}};

    for(int gi=0;ts[gi].g;gi++){
        sprintf(path,"%s/l_s/eigenvalue-3/matrix-test/%s",home,ts[gi].g);
        DIR*dp=opendir(path);if(!dp)continue;
        struct dirent*en;char fp2[1024];
        double*ms=malloc(10000*6*sizeof(double));int n=0;
        while((en=readdir(dp))&&n<10000){
            if(!strstr(en->d_name,".txt"))continue;
            sprintf(fp2,"%s/%s",path,en->d_name);
            FILE*f=fopen(fp2,"r");if(!f)continue;
            double m[6];int ok=1;
            for(int i=0;i<6;i++)if(fscanf(f,"%lf",&m[i])!=1)ok=0;
            fclose(f);if(ok){for(int i=0;i<6;i++)ms[n*6+i]=m[i];n++;}
        }closedir(dp);
        if(!n){free(ms);continue;}

        printf("\n========== %s (n=%d) ==========\n",ts[gi].g,n);

        for(int ei=0;ei<ts[gi].n;ei++){
            double eps=atof(ts[gi].e[ei]);
            int cs[6]={0};double se=0,me=0;int v=0;
            double l[3],lr[3],rr[6];

            for(int i=0;i<n;i++){
                double*m=&ms[i*6];
                int ca=get_case(m,eps);cs[ca]++;
                a2(m,l,eps);
                conv(m,rr);robust(rr,lr);
                double rel=0;
                for(int k=0;k<3;k++){
                    double e=fabs(l[k]-lr[k])/(fabs(lr[k])+1e-300);
                    if(e>rel)rel=e;
                }
                se+=rel;if(rel>me)me=rel;v++;

                sprintf(fp2,"%s/exp_data.csv",base);
                FILE*fw=fopen(fp2,"a");
                if(fw){fprintf(fw,"%s,%s,%d,%.15e\n",ts[gi].g,ts[gi].e[ei],ca,rel);fclose(fw);}
            }

            int red=cs[1]+cs[2]+cs[3]+cs[4];
            printf("eps=%-5s: C1=%3d C2=%3d C3=%3d C4=%3d C5=%3d | red=%3d/%d(%2.0f%%) | avg_err=%.2e max_err=%.2e\n",
                   ts[gi].e[ei],cs[1],cs[2],cs[3],cs[4],cs[5],red,v,100.0*red/v,se/v,me);

            sprintf(fp2,"%s/case_distribution.txt",base);
            fw=fopen(fp2,"a");if(fw){fprintf(fw,"%s,%s,%d,%d,%d,%d,%d,%d,%d,%.0f\n",ts[gi].g,ts[gi].e[ei],cs[1],cs[2],cs[3],cs[4],cs[5],red,v,100.0*red/v);fclose(fw);}
            sprintf(fp2,"%s/error_stats.txt",base);
            fw=fopen(fp2,"a");if(fw){fprintf(fw,"%s,%s,%.2e,%.2e\n",ts[gi].g,ts[gi].e[ei],se/v,me);fclose(fw);}
        }

        int rep=20000;clock_t ca,cb;double ta=0,tc=0;double l[3];
        for(int i=0;i<n;i++){
            double*m=&ms[i*6];
            ca=clock();for(int r=0;r<rep;r++)a2(m,l,1e-10);cb=clock();ta+=(double)(cb-ca);
            ca=clock();for(int r=0;r<rep;r++)a2(m,l,0);cb=clock();tc+=(double)(cb-ca);
        }
        double ta_ns=ta/n/rep/CLOCKS_PER_SEC*1e9,tc_ns=tc/n/rep/CLOCKS_PER_SEC*1e9;
        printf("  Timing: Algo2=%.1fns Cardano=%.1fns ratio=%.2f\n",ta_ns,tc_ns,tc_ns/ta_ns);
        sprintf(fp2,"%s/timing.txt",base);
        fw=fopen(fp2,"a");if(fw){fprintf(fw,"%s,%.1f,%.1f,%.2f\n",ts[gi].g,ta_ns,tc_ns,tc_ns/ta_ns);fclose(fw);}
        free(ms);
    }
    printf("\nDone. Results in ~/l_s/eigenvalue-3/result/\n");
    return 0;
}
