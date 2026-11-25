use std::env;
fn main(){let h:f64=env::args().nth(1).unwrap().parse().unwrap();
let (k,l,w,t,tb,ti)=(180.0f64,0.05f64,0.1f64,0.003f64,200.0f64,40.0f64);
let n=200usize;let p=2.0*(w+t);let a_c=w*t;let m2=h*p/(k*a_c);let dx=l/n as f64;
let mut a=vec![1.0;n+1];let mut b=vec![-(2.0+m2*dx*dx);n+1];let mut cc=vec![1.0;n+1];let mut d=vec![0.0;n+1];let mut x=vec![0.0;n+1];
b[0]=1.0;cc[0]=0.0;d[0]=tb-ti;a[0]=0.0;a[n]=2.0;b[n]=-(2.0+m2*dx*dx+2.0*dx*h/k);cc[n]=0.0;
for i in 1..=n{let f=a[i]/b[i-1];b[i]-=f*cc[i-1];d[i]-=f*d[i-1];}
x[n]=d[n]/b[n];for i in (0..n).rev(){x[i]=(d[i]-cc[i]*x[i+1])/b[i];}
println!("{:.4}",-k*a_c*(-3.0*x[0]+4.0*x[1]-x[2])/(2.0*dx));}