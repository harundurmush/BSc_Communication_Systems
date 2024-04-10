clc;
clear all;
close all;
%%
format long
N=159744;% #bits
SC=64; % #subcarriers
CP=16; % #cycle prefix
SNR=30;
z=zeros(1,SNR+1);
%h=zeros(1,SC);
h=1; %ideal channel
h=[0.04 -0.05 0.07 -0.21 -0.5 0.72 0.36 0 0.21 0.03 0.07 ]; %channel 1
h=[0.407 0.815 0.407]; %channel 2
h=[0.227 0.460 0.688 0.460 0.227]; %cannel 3
smat=zeros(SC,(N/2)/SC);
smatifft=zeros(SC,(N/2)/SC);
xt=zeros(1,SC);
xr=zeros(1,SC);
t=zeros(SC,(N/2)/SC);
xscat=zeros(SC,(N/2)/SC);
xrfftscat=zeros(SC,(N/2)/SC);
rmat=zeros(SC,(N/2)/SC);
xtx=zeros(1,N/2);
ts=zeros(1,N/2);
vectcp=zeros(1,SC+CP);
smatcp=zeros(SC+CP,(N/2)/SC);
vectrcp=zeros(SC,1);
rmatcp=zeros(SC+CP,(N/2)/SC);
I=eye(SC);
W=zeros(1,SC);
%create input chain of bits
x=randi(2,1,N);
x=x-1;
%convert bits into symbols
s=zeros(1,N/2);
b=1;
c=2;
for d=1:N/2
    if x(b)==0 && x(c)==0
        s(d)=1+1j;
    elseif x(b)==0 && x(c)==1
        s(d)=-1+1j;
    elseif x(b)==1 && x(c)==1
        s(d)=-1-1j;
    elseif x(b)==1 && x(c)==0
        s(d)=1-1j;
    end
    b=b+2;
    c=c+2;
end

for d=1:N/2
 smat(d)=s(d);
end
for k=1:(SNR+1)
%ifft
    for e=1:((N/2)/SC)
         xt=smat(:,e);
         xtifft=ifft(xt);
         smatifft(:,e)=xtifft;
    end
    %cyclic prefix
    for e=1:((N/2)/SC)
        smatcp(1:CP,e)=smatifft(SC-CP+1:SC,e);
        smatcp(1+CP:SC+CP,e)=smatifft(1:SC,e);
    end
    %paralel to serial
    for d=1:((N/2)+CP*((N/2)/SC))
        xtx(d)=smatcp(d);
    end
    %variance tx
    potx=0;
    for m=1:((N/2)+CP*((N/2)/SC))
        potx=potx + (abs(xtx(1,m)))^2;
    end
        px=potx/((N/2)+CP*((N/2)/SC));
    %channel
    xtxch=conv(xtx,h);
    trch=zeros(1,(N/2)+CP*((N/2)/SC));
    for d=1:((N/2)+CP*((N/2)/SC))
        trch(d)=xtxch(d);
    end
    %variance after the channel to perform the AWGN
    pot=0;
    for m=1:((N/2)+CP*((N/2)/SC))
     pot=pot + (abs(trch(1,m)))^2;
    end
    ps=pot/((N/2)+CP*((N/2)/SC));
    var=ps/(4*10^((k-1)/10));
    %noise
    nr=randn(1,((N/2)+CP*((N/2)/SC)));
    nr=sqrt(var)*nr;
    ni=randn(1,((N/2)+CP*((N/2)/SC)));
    ni=sqrt(var)*ni;
    ni= 1i*ni;
    n=nr+ni;
    %variance of the noise
    potn=0;
    for m=1:((N/2)+CP*((N/2)/SC))
        potn=potn + (abs(n(1,m))^2);
    end
    pn=potn/((N/2)+CP*((N/2)/SC));
    r = trch + n ;
    %serial to paralel
    for d=1:((N/2)+CP*((N/2)/SC))
        rmatcp(d)=r(d);
    end
    %remove cyclic prefix
    for e=1:((N/2)/SC)
        rmat(:,e)=rmatcp(1+CP:SC+CP,e);
    end
    %coeficients mmse equalizer
    H=fft(h,64);
    Ckmmse=zeros(1,SC);
    for d=1:SC
        Ckmmse(d)= H(d)'/(abs(H(d))^2 + (pn/ps) );
    end
    %coeficients zero forcing equalizer
    H=fft(h,64);
    Ckzf=zeros(1,SC);
    for d=1:SC
        Ckzf(d)= H(d)'/(abs(H(d))^2 + (pn/ps) );
    end
    %fft
    for e=1:((N/2)/SC)
         xr=rmat(:,e);     
         xrfft=fft(xr);
         xrfftscat(:,e)=xrfft;
     for d=1:SC
        xrfft(d)=xrfft(d)*(Ckzf(d)); %equalization
     end
        xscat(:,e)=xrfft;
        t(:,e)=xrfft;
    end
    for d=1:N/2
        ts(d)=t(d);
    end
    f=zeros(1,N);
    fr=real(ts);
    fi=imag(ts);
    %decisor
    for d=1:N/2
         if fr(d)>0
            fr(d)=1;
         elseif fr(d)<0
            fr(d)=-1;
         end
    end
    for d=1:N/2
        if fi(d)>0
            fi(d)=1;
        elseif fi(d)<0
            fi(d)=-1;
        end
    end
    b=1;
    c=2;
    %convert symbols into bits
    for d=1:N/2
         if fr(d)==1 && fi(d)==1
         f(b)=0;
         f(c)=0;
         elseif fr(d)==1 && fi(d)==-1
         f(b)=1;
         f(c)=0;
         elseif fr(d)==-1 && fi(d)==-1
         f(b)=1;
         f(c)=1;
         elseif fr(d)==-1 && fi(d)==1
         f(b)=0;
         f(c)=1;
         end
     b=b+2;
     c=c+2;
end
%BER
d=1;
count=0;
for d=1:N
 if x(d)==f(d)
 else count=count+1;
 end
end
BER=count/N;
z(k)=BER;
end
%theoretical BER
yid=0:0.1:SNR;
z_id=qfunc(sqrt(2*10.^(yid/10)));
%BER plots
y=0:1:SNR;
semilogy(yid,z_id, 'r-')
hold on
semilogy(y,z, 'b*')