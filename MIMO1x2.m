close all;
clear;
N = 1000000;
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
figure
%%
EsN0=0.5*(10.^(snr/10));
Berthe=(1-sqrt(EsN0./(EsN0+1)))-0.25*(1-sqrt(EsN0./(EsN0+1))).^2;
semilogy(snr,Berthe,'gs-','LineWidth',2);
hold on;
%%
semilogy(snr,R_ber,'kp-','LineWidth',2);
axis([0 21 10^-6 1]);
xlabel('SNR(dB)')
ylabel('BER')
legend('BER-theory','BER-1*2');
grid on;
title('单发双收QPSK仿真误码率');
