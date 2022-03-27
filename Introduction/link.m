%Program(1)
clear;
clc;close all;
M=2;   K=4; %M-Modulation,K-channels
N=8;   U=16;  %N-Input symbols per channel,U-Oversampling time
S=randi([0,M-1],1,N*K)-(M-1)/2;%1*NK random integers -(M-1)/2~(M-1)/2
% S=randint(1,N*K,M)-(M-1)/2;
X=reshape(S,K,N);
%% 画离散点信号
figure(1);  %画总的时域信号点
axis([0,N*K+1,0,(K+2)*2*M]);
plot(1:N*K,S+(K+1)*2*M,'.');%N*K symbols in a channel
line([0,N*K+1],[(K+1)*2*M,(K+1)*2*M]);%画线纵坐标为20
hold on;
for kn=1:K %第k个信道，画每个信道的时域信号点
    plot(kn:K:N*K,X(kn,:)+(K-kn+1)*2*M,'.');%N*K symbols in K channels
    line([0,N*K+1],[(K-kn+1)*2*M,(K-kn+1)*2*M]);
end
hold off
title('离散点的分布');legend('总的离散点','第一个信道的离散点','第二个信道的离散点','第三个信道的离散点','第四个信道的离散点');

%% 画傅里叶变换的反变换的时域波形
figure(2);  %画总的傅里叶变换的反变换，相当于是由上面的点变为了连续波波形
axis([0,N*K+1,0,(K+2)*2*M]);
plot(1/U:1/U:N*K,U*real(ifft(fft(S),K*U*N))+(K+1)*2*M);
line([0,N*K+1],[(K+1)*2*M,(K+1)*2*M]);%画线纵坐标为20
hold on;
for kn=1:K  %第k个信道，画每个信道信号傅里叶变换的反变换
    Xkn=X(kn,:);
    Xkt=K*U*real(ifft(fft(Xkn),K*U*N));
    Xkt=[Xkt(K*U*N-(kn-1)*U+1:K*U*N),Xkt(1:K*U*N-(kn-1)*U)];
    plot(1/U:1/U:N*K,Xkt+(K-kn+1)*2*M);hold on;
    line([0,N*K+1],[(K-kn+1)*2*M,(K-kn+1)*2*M]);
end
hold off
title('离散点傅里叶变换再反变换后连续波图像');legend('总的波形图','第一个信道的波形图','第二个信道的波形图','第三个信道的波形图','第四个信道的波形图');

%% 画调制后时域波形
Y=0;
f=0:K-1;
figure(3);  
axis([0,N*K+1,0,(K+2)*2*M]);
plot(1/U:1/U:N*K,U*real(ifft(fft(S),K*U*N))+(K+1)*2*M);%原来的总波形
line([0,N*K+1],[(K+1)*2*M,(K+1)*2*M]);
hold on;
for kn=1:K
    Xkn=X(kn,:);
    Xkt=U*K*real(ifft(fft(Xkn),K*U*N));%和上面的一样
    Xkt=[Xkt(K*U*N-(kn-1)*U+1:K*U*N),Xkt(1:K*U*N-(kn-1)*U)];%和上面的一样
    Ykt=Xkt.*cos(2*pi*f(kn)*(1:K*U*N)/(U)); %乘以载波，调制到通带上（对于第一个通道,f=0,波形何原来一致）
    plot(1/U:1/U:N*K,Ykt+(K-kn+1)*2*M);
    line([0,N*K+1],[(K-kn+1)*2*M,(K-kn+1)*2*M]);
%     Y=Y+Ykt; %这个好像没什么用
end
hold off
%**********************************************************************

%% 计算调制后的BER和SER
%Program(2)
clear;
N=64; % FFT size
i=sqrt(-1);
M=16;  % Modulation level
switch M
    case {2}
        Symbol_Map=[1,-1];
    case {4}
        Symbol_Map=[1+i -1+i 1-i -1-i];
    case {16}
        Symbol_Map=[1+i -1+i 3+i -3+i 1-i -1-i 3-i -3-i 1+3i -1+3i 3+3i -3+3i 1-3i -1-3i 3-3i -3-3i];
    otherwise
        error('Undefined Modulation Level');
end
figure(4);%画16进制的星座图
plot(real(Symbol_Map),imag(Symbol_Map),'*');
axis([-sqrt(M),sqrt(M),-sqrt(M),sqrt(M)]);
Sym_Pow=Symbol_Map*Symbol_Map'; %总功率
Bit_Per_Sym=round(log2(M));     %每个符号多少个比特
Map_Matr=2.^(0:Bit_Per_Sym-1);
Num_Sym=N;          %N个符号
L=16;  % CP length
h_am=[1,0.7,0.5]; % Channel profile，考虑多径效应
h_length=length(h_am);
h=h_am.*(sqrt(0.5)*(randn(1,h_length)+i*randn(1,h_length)));
SNRindB=30;
Noise_Sigma=sqrt(0.5*10^(-SNRindB/10)); %0.5是信道衰减系数，考虑归一化信号功率
Xf_Bit_S=randi([0,1],1,Bit_Per_Sym*Num_Sym);%产生4*64个2进制bit
% Xf_Bit_S=randint(1,Bit_Per_Sym*Num_Sym);
Xf_Bit_P=reshape(Xf_Bit_S,Bit_Per_Sym,Num_Sym);%组合成64个符号，每列是一个符号
SymId0=Map_Matr*Xf_Bit_P+1;     %按每位所代表的值换算成64个10进制数
Xf_Sym=Symbol_Map(SymId0);  %每一个十进制数对应星座图的索引，转变成64个星座点
xt=sqrt(N)*ifft(Xf_Sym,N);      %ifft
xt_CP=[xt(N-L+1:N),xt];         %将上一步傅里叶反变换后xt的后16个数加在xt前面
yt=conv(xt_CP,h);               %和信道进行卷积
Noise=sqrt(Sym_Pow)*Noise_Sigma*(randn(1,N)+i*randn(1,N));%实际每个符号的噪声功率
yk=yt(L+1:L+N)+Noise;   %通过信道输出的总信号
Yf=sqrt(1/N)*fft(yk);   %对应的傅里叶变换
Hf=fft(h,N);            %信道的傅里叶变换
Xd=Yf./Hf;              %输出/信道=原来的信号的（傅里叶变换）?
Xd=Xd(1:Num_Sym);       %这一步好像没什么用
Dec_Matr=abs(ones(M,1)*Xd-Symbol_Map.'*ones(1,Num_Sym));% 64个符号变为16*64
[Distance,SymId]=min(Dec_Matr,[],1);    %取每一列的最小值，及其每一列对应的索引->第几个星座点
Xf_Dec=Symbol_Map(SymId);   %每一个对应的星座点
Xf_Bit_Dec=de2bi(SymId-1,4);%将16进制索引转为4bit形式，每一行为每一个符号
Xf_Bit_Dec=Xf_Bit_Dec';     %每一列为每一个符号
BitError=reshape(Xf_Bit_Dec,1,[])-Xf_Bit_S;     %变为一串2进制比特形式-原来的=比特误差
BER=sum(abs(BitError))/(Num_Sym*Bit_Per_Sym)    %比特误差率
SER=sum(abs(sign(SymId0-SymId)))/Num_Sym        %符号误差率
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%