#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;
int main(int c,char**v){double h=atof(v[1]),k=180,L=0.05,w=0.1,t=0.003,Tb=200,Ti=40;
int N=200;double P=2*(w+t),A=w*t,m2=h*P/(k*A),dx=L/N;
vector<double> a(N+1,1),b(N+1,-(2+m2*dx*dx)),cc(N+1,1),d(N+1,0),x(N+1);
b[0]=1;cc[0]=0;d[0]=Tb-Ti;a[0]=0;a[N]=2;b[N]=-(2+m2*dx*dx+2*dx*h/k);cc[N]=0;
for(int i=1;i<=N;i++){double f=a[i]/b[i-1];b[i]-=f*cc[i-1];d[i]-=f*d[i-1];}
x[N]=d[N]/b[N];for(int i=N-1;i>=0;i--)x[i]=(d[i]-cc[i]*x[i+1])/b[i];
cout.precision(4);cout<<fixed<<-k*A*(-3*x[0]+4*x[1]-x[2])/(2*dx)<<endl;}