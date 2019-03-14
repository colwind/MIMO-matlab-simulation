close all;
clear;
N = 100000;
s10= rand(1,N)>0.5;
s20= rand(1,N)>0.5;
j=sqrt(-1);
%%
ray11=sqrt(0.5)*[randn(1,N/2)+j*randn(1,N/2)];
ray_x11=abs(ray11);
ray12=sqrt(0.5)*[randn(1,N/2)+j*randn(1,N/2)];
ray_x12=abs(ray12);
ray21=sqrt(0.5)*[randn(1,N/2)+j*randn(1,N/2)];
ray_x21=abs(ray21);
ray22=sqrt(0.5)*[randn(1,N/2)+j*randn(1,N/2)];
ray_x22=abs(ray22);
%%
s1= s10*2-1;
s1_i=s1(1:2:N);
s1_q=s1(2:2:N);
s1_QPSK=zeros(1,N/2);
for i=1:(N/2)
    s1_QPSK(i)=s1_i(i)+j*s1_q(i);
end
%%
s2= s20*2-1;
s2_i=s2(1:2:N);
s2_q=s2(2:2:N);
s2_QPSK=zeros(1,N/2);
for i=1:(N/2)
    s2_QPSK(i)=s2_i(i)+j*s2_q(i);
end
%%
y1 = s1_QPSK .* ray_x11+s2_QPSK .* ray_x12;
y2 = s1_QPSK .* ray_x21+s2_QPSK .* ray_x22;

snr = [0 : 20];
EbN0 = 10.^(snr/10);
sp=1;
y_d1=zeros(1,N);
y_d2=zeros(1,N);
for db=1:length(snr)
    np=sp/EbN0(db);%噪声功率
    AWG1=sqrt(0.5*np)*(randn(1,N/2)+j*randn(1,N/2));
    AWG2=sqrt(0.5*np)*(randn(1,N/2)+j*randn(1,N/2));
    y_awg1=y1 + AWG1;
    y_awg2=y2 + AWG2;
    
    for i=1:N/2
        a=[ray_x11(i),ray_x12(i);ray_x21(i),ray_x22(i)];
        b=[y_awg1(i);y_awg2(i)];
        c=a\b;
        y_b1(i)=c(1,1);
        y_b2(i)=c(2,1);
    end    
    
    for k=1:N/2
        y_d1(2*k-1)=real(y_b1(k));
        y_d1(2*k) = imag(y_b1(k));
        y_d2(2*k-1)=real(y_b2(k));
        y_d2(2*k) = imag(y_b2(k));
    end
    y_d1(find(y_d1>0)) = 1;
    y_d1(find(y_d1<=0)) = 0;%解调判决
    y_d2(find(y_d2>0)) = 1;
    y_d2(find(y_d2<=0)) = 0;%解调判决   
    [num_bit1,~]=biterr(s10,y_d1);                     
    R_ber1(db)=num_bit1/N;
    [num_bit2,~]=biterr(s20,y_d2);                     
    R_ber2(db)=num_bit2/N;
end
figure
semilogy(snr,R_ber1,'bx-','LineWidth',2);
hold on;
semilogy(snr,R_ber2,'mx-','LineWidth',2);
hold on;

axis([0 21 10^-3 1]);
xlabel('SNR(dB)')
ylabel('BER')
grid on;
title('双发双收QPSK仿真误码率');

