close all;
clear;
N = 100000;
%%
%建立一个瑞利衰弱信道；采点从0至5，以0.1作为间隔
j=sqrt(-1);
ray=sqrt(0.5)*[randn(1,N/2)+j*randn(1,N/2)];
ray_x=abs(ray);
%ray_x=raylrnd(sqrt(0.5),1,N);
xi=[0:0.1:5];
%x=randn(1,N)+j*randn(1,N);
%计算瑞利衰弱信道x的概率密度函数；
ray_pdf=ksdensity(ray_x,xi);
%作图
figure
plot(xi,ray_pdf,'kp-')
hold on;
grid on;
%%
%%理论值%%
%标准瑞利分布的概率密度函数公式计算
ray_pdf_theory=raylpdf(xi,sqrt(0.5));
%作图
plot(xi,ray_pdf_theory,'gs-');
hold on;
grid on;
legend('瑞利衰弱信道下的pdf','标准瑞利分布的pdf')
xlabel('x');
ylabel('fx');
title('瑞利衰弱信道与标准瑞利分布的pdf曲线比较')
%%
%单发单收，QPSK调制
%星座点映射
s10= rand(1,N)>0.5;
s1= s10*2-1;
s1_i=s1(1:2:N);
s1_q=s1(2:2:N);
s1_QPSK=zeros(1,N/2);
for i=1:(N/2)
    s1_QPSK(i)=s1_i(i)+j*s1_q(i);
end
%scatterplot(s1_QPSK);
%加信道
s1_qpsk_ray = s1_QPSK * ray_x;%rayleigh
%scatterplot(s1_qpsk_ray);

snr = [0 : 20];
EbN0 = 10.^(snr/10);
sp=1;
y_dqpsk=zeros(1,N);
for db=1:length(snr)
    np=sp/EbN0(db);%噪声功率
    AWG=sqrt(0.5*np)*(randn(1,N/2)+j*randn(1,N/2));
    s1_qpsk_ray_awg=s1_qpsk_ray + AWG;
    y = s1_qpsk_ray_awg;
    y_b = y / ray_x;
    
    for k=1:N/2
        y_dqpsk(2*k-1)=real(y_b(k));
        y_dqpsk(2*k)=imag(y_b(k));
    end
    y_dqpsk(find(y_dqpsk>=0)) = 1;
    y_dqpsk(find(y_dqpsk<0)) = 0;%解调判决
    
 %   R_ber(db)=(sum(abs((y_dqpsk-s1))))/N;
  %  R_ber(db)=size(find([y_dqpsk-s1],2));    
   [num_bit,~]=biterr(s10,y_dqpsk);                     
   R_ber(db)=num_bit/N; 
end
EsN0=0.5*(10.^(snr/10));
Berthe=(1-sqrt(EsN0./(EsN0+1)))-0.25*(1-sqrt(EsN0./(EsN0+1))).^2;

figure
semilogy(snr,R_ber,'bx-','LineWidth',2);
hold on;
semilogy(snr,Berthe,'gs-','LineWidth',2);
axis([0 21 10^-3 1]);

hold on;
grid on;














