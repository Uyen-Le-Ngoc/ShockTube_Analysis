// icc code.cpp -o code $MKLP
#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <stdlib.h>
#include <string>
#include "omp.h"
#include <mkl.h>
using namespace std;

const int N=1000; //step of time
const int Nx=10000; //step of x
double L=4; // length of tube
double dt=1e-3, dx=L/Nx;
//initial state
double ga=1.4;
double r01=0.125,r05=1;
double p01=0.1,p05=1;
double G=(ga-1)/(ga+1);
double alpha=2*ga/(ga-1);
double x0=0; //position of interface

double p[Nx], p1,p2,p3,p4,p5; //pressure
double r[Nx], r1,r2,r3,r4,r5; // Density rpo of each area
double x[Nx];
double u[N], u1=0,u2,u3,u4,u5=0; //velocity
double vs,vc,vr0,vr1; //interface wave velocity
double P,t; //Finding by itetation
double m1,m2,m3,m4; // interface position
double c[Nx], c1=sqrt(ga*p01/r01),c2,c3,c4,c5=sqrt(ga*p05/r05);

void Writedata(int n)
{
    char data[256];
    sprintf(data, "out_%d", n);
    ofstream ofs(data);
    for (int i=0; i < Nx; i++){
        x[i]=x[0]+i*dx;
      //ofs << "x"<< ' ' << "density" << ' ' << "velocity" << "pressure" << endl; 
        ofs << x[i]  << ' ' << r[i] << ' ' << u[i] << ' ' << p[i] <<  endl;
    }
        ofs << endl;
}

int main(){
//value of P
for (int i=0; i<5000; i++)
{
    P=0.001*i;
    double Z=sqrt(2/(ga*((ga-1)+(ga+1)*P)))*(P-1)-2*c1/(c5*(ga-1))*(pow((p01/p05),1/alpha)-pow(P,1/alpha))+(u1-u5)/c5;
    double Z1=c1*sqrt(2/ga)*(P-1)/sqrt(ga-1+(ga+1)*P)-2*c5*(1-pow((P*p01/p05),(1/alpha)))/(ga-1);
//    cout << P << ' ' << Z << ' ' << Z1 << endl;
    if (abs(Z1)<1e-3) {    
    break;
    }
}
x[0]=-L/2;

for (int n=0; n<=N; n++)
{
if (n%10==0){Writedata(n/10);}
for (int i=0; i<Nx; i++)
{
    x[i]=x[0]+i*dx;
    vs=u1+c1*sqrt((ga-1+(ga+1)*P)/(2*ga));
    m1=x0+vs*t;
    //REGION1------x>x0+vs*t
    if (L/2>x[i] && x[i]>x0+vs*t)
    {
    r[i]=r01;
    p[i]=p01;
    u[i]=u1;
    c[i]=sqrt(ga*p[i]/r[i]);
    }

    //REGION5------x<x0+(u5-c5)*t
    u5=0;
    m4=x0+(u5-c5)*t;
    if (-L/2<x[i] && x[i]<x0+(u5-c5)*t){
    r[i]=r05;
    p[i]=p05;
    u[i]=0;
    }
    
    vc=u5+2*c5*(1-pow((P*p01/p05),(1/alpha)))/(ga-1);
    //REGION2------x0+vc*t<x<=x0+vs*t
    m2=x0+vc*t;
    if (x[i]>x0+vc*t && x[i]<=x0+vs*t){
    r[i]=r01*(P+G)/(1+G*P);
    p[i]=P*p01;
    u[i]=vs+c1*sqrt(2/ga)*(P-1)/sqrt(ga-1+(ga+1)*P);
    }
    c2=c1*sqrt(P*(1+G*P)/(G+P));

    //REGION3------x0+((ga+1)*vc/2-c5-(ga-1)*u5/2)*t<x<=x0+vc*t
    m3=x0+((ga+1)*vc/2-c5-(ga-1)*u5/2)*t;
    if (x[i]>x0+((ga+1)*vc/2-c5-(ga-1)*u5/2)*t && x[i]<=x0+vc*t){
    p[i]=P*p01;
    u[i]=vc;
    r[i]=r05*(pow(P*p01/p05,1/ga));
    }

    //REGION4------x0+(u5-c5)*t<=x<x0+((ga+1)*vc/2-c5-(ga-1)*u5/2)*t
    if (x0+(u5-c5)*t <= x[i] && x[i]<x0+((ga+1)*vc/2-c5-(ga-1)*u5/2)*t)
    {
    u[i]=2*(c5+(x[i]-x0)/t+(ga-1)*u5/2)/(ga+1);
    c4=c5-(ga-1)*(u[i]-u5)/2;
    p[i]=p05*pow(c4/c5,alpha);
    r[i]=r05*pow(p[i]/p05,1/ga);
    double vr=-c5;
    }
//    cout << m1 << ' ' << m2 << ' '<< m3 << ' ' << m4 <<endl;
}
    t=n*dt;
}
}