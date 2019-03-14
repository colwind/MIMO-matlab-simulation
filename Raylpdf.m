clc;
clear;
N=50000;
M=randn(1,N)+j*randn(1,N);
r=abs(M);
[f,xi]=ksdensity(r);
subplot(1,2,1)
plot(xi,f);
legend('瑞利衰落pdf曲线');
hold on
x=-2:0.1:5;
p=raylpdf(x,1);
subplot(1,2,2)
plot(x,p,'r-');
legend('标准瑞利衰落pdf');


    