close all;
clear;
N = 1000000;
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
s1_qpsk_ray = s1_QPSK .* ray_x;%rayleigh
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
    y_b = y./ ray_x;
    
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

figure
semilogy(snr,R_ber,'bx-','LineWidth',2);
hold on;

s10= rand(1,N)>0.5;
j=sqrt(-1);
%%
ray11=sqrt(0.5)*[randn(1,N/2)+j*randn(1,N/2)];
ray_x11=abs(ray11);
ray21=sqrt(0.5)*[randn(1,N/2)+j*randn(1,N/2)];
ray_x21=abs(ray21);
%%
s1= s10*2-1;
s1_i=s1(1:2:N);
s1_q=s1(2:2:N);
s1_QPSK=zeros(1,N/2);
for i=1:(N/2)
    s1_QPSK(i)=s1_i(i)+j*s1_q(i);
end

%%
y1 = s1_QPSK .* ray_x11;
y2 = s1_QPSK .* ray_x21;

snr = [0 : 20];
EbN0 = 10.^(snr/10);
sp=1;
y_d=zeros(1,N);
for db=1:length(snr)
    np=sp/EbN0(db);%噪声功率
    AWG1=sqrt(0.5*np)*(randn(1,N/2)+j*randn(1,N/2));
    AWG2=sqrt(0.5*np)*(randn(1,N/2)+j*randn(1,N/2));
    
    y_awg1=y1 + AWG1;
    y_awg2=y2 + AWG2;
    
    y_b= conj(ray_x11) .* y_awg1 + conj(ray_x21) .* y_awg2 ;
    
    for k=1:N/2
        y_d(2*k-1)=real(y_b(k));
        y_d(2*k)=imag(y_b(k));
    end
    y_d(find(y_d>=0)) = 1;
    y_d(find(y_d<0)) = 0;%解调判决
    [num_bit,~]=biterr(s10,y_d);                     
    R_ber(db)=num_bit/N;
end
%%
semilogy(snr,R_ber,'kp-','LineWidth',2);
N = 5000000;
s10= rand(1,N)>0.5;
%%
ray1=sqrt(0.5)*[randn(1,N/2)+j*randn(1,N/2)];
ray_x1=abs(ray1);
ray2=sqrt(0.5)*[randn(1,N/2)+j*randn(1,N/2)];
ray_x2=abs(ray2);
ray3=sqrt(0.5)*[randn(1,N/2)+j*randn(1,N/2)];
ray_x3=abs(ray3);
%%
s1= s10*2-1;
s1_i=s1(1:2:N);
s1_q=s1(2:2:N);
s1_QPSK=zeros(1,N/2);
for i=1:(N/2)
    s1_QPSK(i)=s1_i(i)+j*s1_q(i);
end

%%
y1 = s1_QPSK .* ray_x1;
y2 = s1_QPSK .* ray_x2;
y3 = s1_QPSK .* ray_x3;
%%
snr = [0 : 20];
EbN0 = 10.^(snr/10);
sp=1;
y_d=zeros(1,N);
for db=1:length(snr)
    np=sp/EbN0(db);%噪声功率
    AWG1=sqrt(0.5*np)*(randn(1,N/2)+j*randn(1,N/2));
    AWG2=sqrt(0.5*np)*(randn(1,N/2)+j*randn(1,N/2));
    AWG3=sqrt(0.5*np)*(randn(1,N/2)+j*randn(1,N/2));
    y_awg1=y1 + AWG1;
    y_awg2=y2 + AWG2;
    y_awg3=y3 + AWG3;
    y_b= conj(ray_x1) .* y_awg1 + conj(ray_x2) .* y_awg2 + conj(ray_x3) .* y_awg3 ;
    
    for k=1:N/2
        y_d(2*k-1)=real(y_b(k));
        y_d(2*k)=imag(y_b(k));
    end
    y_d(find(y_d>=0)) = 1;
    y_d(find(y_d<0)) = 0;%解调判决
    [num_bit,~]=biterr(s10,y_d);                     
    R_ber(db)=num_bit/N;
end
%%
semilogy(snr,R_ber,'mx-','LineWidth',2);
axis([0 21 10^-7 0.5]);
xlabel('SNR(dB)')
ylabel('BER')
legend('BER-1*1','BER-1*2','BER-1*3');
grid on;
title('误码率比较');









