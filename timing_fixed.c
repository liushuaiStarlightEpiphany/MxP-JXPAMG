#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dirent.h>
#include <time.h>

#define PI 3.141592653589793

static int sc(double x){return(x>0)-(x<0);}
static int cd(const void*a,const void*b){double d=*(const double*)a-*(const double*)b;return(d>0)-(d<0);}

static void a2(const double m[6],double l[3],double eps){
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

static void jacobi(const double m[6],double l[3]){
    double A[9]={m[0],m[3],m[4],m[3],m[1],m[5],m[4],m[5],m[2]};
    for(int it=0;it<50;it++){
        double mx=0;int p=0,q=1;
        for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)if(fabs(A[i*3+j])>mx){mx=fabs(A[i*3+j]);p=i;q=j;}
        if(mx<1e-30)break;double th=(A[q*3+q]-A[p*3+p])/(2*A[p*3+q]);
        double t=sc(th)/(fabs(th)+sqrt(1+th*th)),c=1/sqrt(1+t*t),s=c*t,tau=s/(1+c);
        double ap=A[p*3+p],aq=A[q*3+q],apq=A[p*3+q];
        A[p*3+p]=ap-t*apq;A[q*3+q]=aq+t*apq;A[p*3+q]=A[q*3+p]=0;
        for(int r=0;r<3;r++)if(r!=p&&r!=q){double arp=A[r*3+p],arq=A[r*3+q];A[r*3+p]=A[p*3+r]=arp-s*(arq+tau*arp);A[r*3+q]=A[q*3+r]=arq+s*(arp-tau*arq);}
    }l[0]=A[0];l[1]=A[4];l[2]=A[8];qsort(l,3,sizeof(double),cd);
}

static void qrf(const double m[6],double l[3]){
    double A[9]={m[0],m[3],m[4],m[3],m[1],m[5],m[4],m[5],m[2]};
    for(int it=0;it<50;it++){
        double sh=A[8];for(int i=0;i<3;i++)A[i*3+i]-=sh;
        for(int i=0;i<2;i++){double a=A[i*3+i],b=A[(i+1)*3+i],r=hypot(a,b),c=a/r,s=-b/r;
            for(int j=i;j<3;j++){double x=A[i*3+j],y=A[(i+1)*3+j];A[i*3+j]=c*x-s*y;A[(i+1)*3+j]=s*x+c*y;}
            for(int j=0;j<3;j++){double x=A[j*3+i],y=A[j*3+i+1];A[j*3+i]=c*x-s*y;A[j*3+i+1]=s*x+c*y;}
        }for(int i=0;i<3;i++)A[i*3+i]+=sh;
    }l[0]=A[0];l[1]=A[4];l[2]=A[8];qsort(l,3,sizeof(double),cd);
}

typedef struct{const char*g;double eps;int n;}TG;

int main(){
    char home[256];sprintf(home,"%s",getenv("HOME"));
    char p[1024];

    TG ts[]={
        {"diag",1e-4,500},{"ps1",1e-4,500},{"ps2",1e-4,500},{"ps3",1e-4,500},
        {"general",1e-4,500},{"dual",1e-4,300},
        {"sweep/2.5e-3",1e-2,100},{"sweep/2.5e-4",1e-3,100},
        {"sweep/2.5e-5",1e-4,100},{"sweep/2.5e-6",1e-5,100},{"sweep/2.5e-7",1e-6,100},
        {NULL,0,0}
    };

    printf("%-16s %8s %8s %8s %8s  (%s)\n","group","Algo2","Cardano","Jacobi","QR","eps for a2");
    for(int gi=0;ts[gi].g;gi++){
        sprintf(p,"%s/l_s/eigenvalue-3/matrix-test/%s",home,ts[gi].g);
        DIR*dp=opendir(p);if(!dp)continue;struct dirent*en;char fp2[1024];
        double*ms=malloc(ts[gi].n*6*sizeof(double));int n=0;
        while((en=readdir(dp))&&n<ts[gi].n){
            if(!strstr(en->d_name,".txt"))continue;
            sprintf(fp2,"%s/%s",p,en->d_name);FILE*f2=fopen(fp2,"r");if(!f2)continue;
            double m[6];int ok=1;
            for(int i=0;i<6;i++)if(fscanf(f2,"%lf",&m[i])!=1)ok=0;
            fclose(f2);if(ok){for(int i=0;i<6;i++)ms[n*6+i]=m[i];n++;}
        }closedir(dp);
        if(!n){free(ms);continue;}

        double eps=ts[gi].eps;
        int rep=50000;clock_t ca,cb;double l[3];
        double ta=0,tc=0,tj=0,tq=0;

        for(int i=0;i<n;i++){
            double*m=&ms[i*6];
            ca=clock();for(int r=0;r<rep;r++)a2(m,l,eps);cb=clock();ta+=(double)(cb-ca);
            ca=clock();for(int r=0;r<rep;r++)a2(m,l,0);cb=clock();tc+=(double)(cb-ca);
            ca=clock();for(int r=0;r<rep/4;r++)jacobi(m,l);cb=clock();tj+=(double)(cb-ca);
            ca=clock();for(int r=0;r<rep/4;r++)qrf(m,l);cb=clock();tq+=(double)(cb-ca);
        }
        double t_a2=ta/n/rep/CLOCKS_PER_SEC*1e9;
        double t_ca=tc/n/rep/CLOCKS_PER_SEC*1e9;
        double t_ja=tj/n/(rep/4)/CLOCKS_PER_SEC*1e9;
        double t_qr=tq/n/(rep/4)/CLOCKS_PER_SEC*1e9;
        printf("%-16s %7.1f %7.1f %7.1f %7.1f  (eps=%.0e)\n",ts[gi].g,t_a2,t_ca,t_ja,t_qr,eps);
        free(ms);
    }
    return 0;
}
