close all  
clc 
clear all  
N=100000; 
xn0 = randi([0 1],1,N);%原序列
xn1 = randi([0 1],1,N);

for SNR=0:20
    
for k=1:2:N  %QPSK调制，将xn每两个分成一个虚数，mn为虚数列
    t=(k+1)/2;
    if xn0(k)==0 && xn0(k+1)==0
        mn0(t)=-sqrt(2)/2-(sqrt(2)/2)*j;
    elseif xn0(k)==0 && xn0(k+1)==1
        mn0(t)=-sqrt(2)/2+(sqrt(2)/2)*j;
    elseif xn0(k)==1 && xn0(k+1)==0
        mn0(t)=sqrt(2)/2-(sqrt(2)/2)*j;
    elseif xn0(k)==1 && xn0(k+1)==1
        mn0(t)=sqrt(2)/2+(sqrt(2)/2)*j;
    end
end
for k=1:2:N  %QPSK调制，将xn每两个分成一个虚数，mn为虚数列
    t=(k+1)/2;
    if xn1(k)==0 && xn1(k+1)==0
        mn1(t)=-sqrt(2)/2-(sqrt(2)/2)*j;
    elseif xn1(k)==0 && xn1(k+1)==1
        mn1(t)=-sqrt(2)/2+(sqrt(2)/2)*j;
    elseif xn1(k)==1 && xn1(k+1)==0
        mn1(t)=sqrt(2)/2-(sqrt(2)/2)*j;
    elseif xn1(k)==1 && xn1(k+1)==1
        mn1(t)=sqrt(2)/2+(sqrt(2)/2)*j;
    end
end

h00=Rayleigh(N/2);
h01=Rayleigh(N/2);
h10=Rayleigh(N/2);
h11=Rayleigh(N/2);

snr=1/(10^(SNR/10)); %snr=2v^2
noise1=sqrt(snr/4)*(randn(1,N/2));
noise2=sqrt(snr/4)*(randn(1,N/2))*j;
noise3=sqrt(snr/4)*(randn(1,N/2));
noise4=sqrt(snr/4)*(randn(1,N/2))*j;

y1=mn0.*h00+mn1.*h01+noise1+noise2;
y2=mn0.*h10+mn1.*h11+noise3+noise4;

for i=1:N/2
    a=[h00(i),h01(i);h10(i),h11(i)];
    b=[y1(i);y2(i)];
    c=a\b;
    n1(i)=c(1,1);
    n2(i)=c(2,1);
end

for i=1:2:N   %QPSK解调
    t=(i+1)/2;
    r1(i)=real(n1(t));
end
for i=2:2:N    
    t=i/2;
    r1(i)=imag(n1(t));
end

for i=1:2:N   %QPSK解调
    t=(i+1)/2;
    r2(i)=real(n2(t));
end
for i=2:2:N    
    t=i/2;
    r2(i)=imag(n2(t));
end

for k=1:N 
    if r1(k)>0
        Q1(k) = 1;    %Q为经过调制解调收到的序列
    else 
        Q1(k) = 0;
    end
end
for k=1:N 
    if r2(k)>0
        Q2(k) = 1;    %Q为经过调制解调收到的序列
    else 
        Q2(k) = 0;
    end
end

pbit1(SNR+1)=(sum(abs((xn0-Q1))))/N;%QPSK误比特率计算 
pbit2(SNR+1)=(sum(abs((xn1-Q2))))/N;%QPSK误比特率计算 
end
r=0:20;
semilogy(r,pbit1,'b-s',r,pbit2,'r-*')
xlabel('SNR(dB)')
ylabel('Pe')
grid on;
legend('2发2收QPSK仿真误码率');
    