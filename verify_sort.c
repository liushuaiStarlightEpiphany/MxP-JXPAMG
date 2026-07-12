#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
static int cd(const void*a,const void*b){double d=*(const double*)a-*(const double*)b;return(d>0)-(d<0);}
void robust(const double A[6],double lambda[3]);
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
    else{double phi=acos(dq/sqrt(-dp*dp*dp))/3,r=2*sqrt(-dp);l[0]=r*cos(phi)+t3;l[1]=r*cos(phi+2*M_PI/3)+t3;l[2]=r*cos(phi+4*M_PI/3)+t3;}
    qsort(l,3,sizeof(double),cd);
}
int main(){
    double m[6];FILE*fp=fopen("/home/export/base/sc100286/sc100286/l_s/eigenvalue-3/matrix-test/general/matrix_00000.txt","r");
    for(int i=0;i<6;i++)fscanf(fp,"%lf",&m[i]);fclose(fp);
    double l[3],lr[3],rr[6];
    a2(m,l,1e-4);
    rr[0]=m[0];rr[1]=m[3];rr[2]=m[4];rr[3]=m[1];rr[4]=m[5];rr[5]=m[2];
    robust(rr,lr);qsort(lr,3,sizeof(double),cd);
    printf("m = [%.2e %.2e %.2e | %.2e %.2e %.2e]\n",m[0],m[1],m[2],m[3],m[4],m[5]);
    printf("algo2: %.6e %.6e %.6e\n",l[0],l[1],l[2]);
    printf("robust: %.6e %.6e %.6e\n",lr[0],lr[1],lr[2]);
    double e0=fabs(l[0]-lr[0])/fabs(lr[0]+1e-300);
    double e1=fabs(l[1]-lr[1])/fabs(lr[1]+1e-300);
    double e2=fabs(l[2]-lr[2])/fabs(lr[2]+1e-300);
    printf("errs: %.6e %.6e %.6e\n",e0,e1,e2);
    int ca=0;
    double m0=fmin(fmin(fabs(m[0]),fabs(m[1])),fabs(m[2]));
    double r1=2*fabs(m[3])/m0,r2=2*fabs(m[4])/m0,r3=2*fabs(m[5])/m0;
    if(r1+r2+r3<1e-4)ca=1;else if(r1+r2<1e-4)ca=2;else if(r1+r3<1e-4)ca=3;else if(r2+r3<1e-4)ca=4;else ca=5;
    printf("case=%d r1=%.2e r2=%.2e r3=%.2e r12+r3=%.2e\n",ca,r1,r2,r3,r1+r2+r3);
    return 0;
}
