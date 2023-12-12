scf(0);clf(0);clear;clc;   //Prac Q1
h=1973;m=0.511e6;e=3.795
rmin=1e-10;rmax=20;
n=1000;
r=linspace(rmin,rmax,n);
d=r(2)-r(1);
t=-(h^2)/(2*m*(d^2));
V=zeros(n,n);
K=-2*eye(n,n)+diag(ones(n-1,1),1) +diag(ones(n-1,1),-1);
for i=1:n
    V(i,i)=-(e^2)/r(i);
end
H=t*K+V;
[U1,EV]=spec(H);
E=diag(EV);
disp("Ground state energy : "+string(E(2))+" eV")
disp("1st Excited state energy : "+string(E(3))+" eV")
disp("2nd Excited state energy : "+string(E(4))+" eV")
disp("3rd Excited state energy : "+string(E(5))+" eV")
plot(r', [U1(:,2),U1(:,3),U1(:,4),U1(:,5)],'linewidth',2)
legend('Ground state','First Excited state','Second excited state','Third excited state',4)
title("Plot of wavefunction for e^2/r Potential",'fontsize',3)
xlabel("r(A)",'fontsize',2)
ylabel("Wavefunction",'fontsize',2)
xgrid()
a=gca();
a.x_location="origin";
a.y_location="origin";
