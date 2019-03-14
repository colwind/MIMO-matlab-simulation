close all;
clear all;

N=10^5;

%%  ������һ������źŲ�����qpsk����
random_signal=rand(1,N)>0.5;
r_s=zeros(1,N);
r_s=2*random_signal-1;
rs1=r_s(1:2:N);
rs2=r_s(2:2:N);
rs_qpsk=rs1+j*rs2;

%%  �����ڶ�������źŲ�����qpsk����
random_signal_1=rand(1,N)>0.5;
r_s_1=zeros(1,N);
r_s_1=2*random_signal_1-1;
rs1_1=r_s_1(1:2:N);
rs2_1=r_s_1(2:2:N);
rs_qpsk_1=rs1_1+j*rs2_1;
%%  rayleigh channel
j=sqrt(-1);
rayleigh_channel=1/sqrt(2)*(randn(1,N/2)+j*randn(1,N/2));
rayleigh_channel1=1/sqrt(2)*(randn(1,N/2)+j*randn(1,N/2));
rayleigh_channel2=1/sqrt(2)*(randn(1,N/2)+j*randn(1,N/2));
rayleigh_channel3=1/sqrt(2)*(randn(1,N/2)+j*randn(1,N/2));


%%  ����˥�䣬��˹����
snr=0:20;
SP=sqrt(2);
Eb_N0=10.^(snr/10);

rayleigh_signal=rs_qpsk.*rayleigh_channel+rs_qpsk_1.*rayleigh_channel1;
rayleigh_signal1=rs_qpsk.*rayleigh_channel2+rs_qpsk_1.*rayleigh_channel3;

    %ysiso1=rs_qpsk.*rayleigh_channel;
for db=1:length(snr)
    NP=SP/Eb_N0(db);    %��������
    WGN=sqrt(0.5*NP)*(randn(1,N/2)+j*randn(1,N/2));    %��˹������
    WGN1=sqrt(0.5*NP)*(randn(1,N/2)+j*randn(1,N/2));
    y=rayleigh_signal+WGN;
    y1=rayleigh_signal1+WGN1;
%          ysiso=(ysiso1+WGN)./rayleigh_channel;
    for i=1:N/2
        a=[rayleigh_channel(i),rayleigh_channel1(i);rayleigh_channel2(i),rayleigh_channel3(i)];
        b=[y(i);y1(i)];
        c=a\b;
        rec(i)=c(1,1);
        rec1(i)=c(2,1);
    end
    dqpsk=zeros(1,N);
    dqpsk_1=zeros(1,N);
          % dqpsksiso=zeros(1,N);
    %%  ������о�
    for k=1:N/2;
        dqpsk(2*k-1)=real(rec(k));
        dqpsk(2*k)=imag(rec(k));
        dqpsk_1(2*k-1)=real(rec1(k));
        dqpsk_1(2*k)=imag(rec1(k)); 
%            dqpsksiso(2*k-1)=real(ysiso(k));
%           dqpsksiso(2*k)=imag(ysiso(k));   
    end
    dqpsk(find(dqpsk>0))=1;
    dqpsk(find(dqpsk<=0))=0;
    dqpsk_1(find(dqpsk_1>0))=1;
    dqpsk_1(find(dqpsk_1<=0))=0;
         %dqpsksiso(find(dqpsksiso>0))=1;
        % dqpsksiso(find(dqpsksiso<=0))=0;
    
    %%  ͳ��������
    E(db)=size(find([dqpsk-random_signal]),2);
    E1(db)=size(find([dqpsk_1-random_signal_1]),2);
         % Esiso(db)=size(find([dqpsksiso-random_signal]),2);
end
%% ��ͼ
RBER=E/N;
RBER1=E1/N;
%     RBERsiso=Esiso/N; 
EsN0=0.5*(10.^(snr/10));
RBER2=(1-sqrt(EsN0./(EsN0+1)))-0.25*(1-sqrt(EsN0./(EsN0+1))).^2;

figure

semilogy(snr,RBER,'bx-','LineWidth',2);
hold on;
semilogy(snr,RBER1,'mx-','LineWidth',2);
hold on;

xlabel('SNR(dB)')
ylabel('BER')
grid on;
title('˫��˫��QPSK����������');

axis([0 20 10^-3 0.5] );
grid on;

%legend('QPSK-Simulation1','QPSK-Simulation2','QPSK-Theory');
legend('Simulation1','Simulation2');
