package main
import ("fmt";"os";"strconv")
func main(){h,_:=strconv.ParseFloat(os.Args[1],64)
k,L,w,t,Tb,Ti:=180.0,0.05,0.1,0.003,200.0,40.0
N:=200;P:=2*(w+t);A:=w*t;m2:=h*P/(k*A);dx:=L/float64(N)
a:=make([]float64,N+1);b:=make([]float64,N+1);cc:=make([]float64,N+1);d:=make([]float64,N+1);x:=make([]float64,N+1)
for i:=range a{a[i]=1;b[i]=-(2+m2*dx*dx);cc[i]=1}
b[0]=1;cc[0]=0;d[0]=Tb-Ti;a[0]=0;a[N]=2;b[N]=-(2+m2*dx*dx+2*dx*h/k);cc[N]=0
for i:=1;i<=N;i++{f:=a[i]/b[i-1];b[i]-=f*cc[i-1];d[i]-=f*d[i-1]}
x[N]=d[N]/b[N];for i:=N-1;i>=0;i--{x[i]=(d[i]-cc[i]*x[i+1])/b[i]}
fmt.Printf("%.4f\n",-k*A*(-3*x[0]+4*x[1]-x[2])/(2*dx))}