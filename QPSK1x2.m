close all  
clc 
clear all  

N=300000; 
xn = randi([0 1],1,N);%ԭ����

for SNR=0:20
    
for k=1:2:N  %QPSK���ƣ���xnÿ�����ֳ�һ��������mnΪ������
    t=(k+1)/2;
    if xn(k)==0 && xn(k+1)==0
        mn(t)=-sqrt(2)/2-(sqrt(2)/2)*j;
    elseif xn(k)==0 && xn(k+1)==1
        mn(t)=-sqrt(2)/2+(sqrt(2)/2)*j;
    elseif xn(k)==1 && xn(k+1)==0
        mn(t)=sqrt(2)/2-(sqrt(2)/2)*j;
    elseif xn(k)==1 && xn(k+1)==1
        mn(t)=sqrt(2)/2+(sqrt(2)/2)*j;
    end
end

h0=Rayleigh(N/2);
sn=mn.*h0;
h1=Rayleigh(N/2);
pn=mn.*h1;


snr=1/(10^(SNR/10)); %snr=2v^2
noise1=sqrt(snr/4)*(randn(1,N/2));
noise2=sqrt(snr/4)*(randn(1,N/2))*j;
sn=sn+noise1+noise2;%�����˹������
noise3=sqrt(snr/4)*(randn(1,N/2));
noise4=sqrt(snr/4)*(randn(1,N/2))*j;
pn=pn+noise3+noise4;%�����˹������


y=conj(h0).*sn+conj(h1).*pn;

%sn=sn.*(conj(Ray)./(Ray.*conj(Ray)));

for i=1:2:N   %QPSK���
    t=(i+1)/2;
    r(i)=real(y(t));
end
for i=2:2:N    
    t=i/2;
    r(i)=imag(y(t));
end

for k=1:N 
    if r(k)>0
        Q(k) = 1;    %QΪ�������ƽ���յ�������
    else 
        Q(k) = 0;
    end
end

pbit(SNR+1)=(sum(abs((xn-Q))))/N;%QPSK������ʼ��� 
     
end

r=0:20
semilogy(r,pbit,'b-s')
xlabel('SNR(dB)')
ylabel('Pe')
grid on;
legend('QPSK����������');
    