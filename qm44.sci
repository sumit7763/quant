clc;clear;clf();
h=197.3;m=940e6;D=0.755501;
rmin=1e-3;rmax=2
n=200;a=1.44;r0=0.131349;mu=0.99946*m
r=linspace(rmin,rmax,n);
d=r(2)-r(1);
t=-(h^2)/(2*mu*(d^2));
V=zeros(n,n);
U=-2*eye(n,n)+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
for i=1:n
    rd=(r(i)-r0)/r(i);
    V(i,i)=D*(exp(-2*a*rd)-exp(-a*rd));
end
H=t*U+V;
[U1,EV]=spec(H);
E=diag(EV);
disp("Ground state energy : "+string(E(1))+"eV")
plot(r',U1(:,1),'linewidth',2)
title("Plot of wavefunction for V=D(exp(-2ar)-exp(-ar))")
xlabel("r(A)--------->")
ylabel("psi----------->")
xgrid()
legend("Ground State","First excited state",-4)
