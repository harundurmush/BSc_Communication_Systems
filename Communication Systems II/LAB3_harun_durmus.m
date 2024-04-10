clc;
clear all;
close all;
%% step a
Fs = 3000;
Td = 0.4;
Ts = 1/Fs;
t=0:Ts:(Td-Ts);
%% step b
N = 4;
Tb=Td/N;
%% step c
s0=[zeros(1,Tb/(4*Ts)) ones(1,Tb/(4*Ts)) ones(1,Tb/(4*Ts)) zeros(1,Tb/(4*Ts))];
s1=[(-1)*(ones(1,Tb/(3*Ts))) zeros(1,Tb/(3*Ts)) ones(1,Tb/(3*Ts))];

figure(1)
subplot(211)
plot(0:Ts:Tb-Ts,s0,"color","#A2142F");
xlabel("time (s)");
ylabel("s0(t) (V)");
legend("s0(t)");
title("s_0(t)");
subplot(212)
plot(0:Ts:Tb-Ts,s1,"color","#EDB120");
xlabel("time (s)");
ylabel("s1(t) (V)");
legend("s1(t)");
title("s_1(t)");
%% step d
b=[1 0 1 0];
st = [];
for k=1:length(b)
    if (b(k)==0)
        st =[st s0];
    else
        st =[st s1];
    end
end

figure(2)
plot(t,st,"color","#D95319");
xlabel("time (s)");
ylabel("s(t) (V)");
legend("s(t)");
title("Transmitted signal s(t)");
%% step e
Pavg= sum(abs(st).^2)/length(st);
display(Pavg);
%% step f
% r(t) = s(t) + n(t)
%% step g
snr_db=[15 0];
snr_lin(1)=10^(0.1*snr_db(1));
snr_lin(2)=10^(0.1*snr_db(2));
var(1)=Pavg/snr_lin(1);
var(2)=Pavg/snr_lin(2);
awgn1=sqrt(var(1))*randn(1,length(st));
awgn2=sqrt(var(2))*randn(1,length(st));
rt1 = st + awgn1;
rt2 = st + awgn2;

figure(3)
plot(t,rt1,"color","#0072BD");
xlabel("time (s)");
ylabel("r1(t) (V)");
legend("r1(t)");
title("Received signal over 15 dB SNR");

figure(4)
plot(t,rt2,"color","#A2142F");
xlabel("time (s)");
ylabel("r2(t) (V)");
legend("r2(t)");
title("Received signal over 0 dB SNR");
%% step h
Wb=Tb/Ts; % step i
r10=[];
r11=[];
r20=[];
r21=[];
for k=1:N
    sum10=0;
    sum11=0;
    sum20=0;
    sum21=0;
    for n=((k-1)*Wb+1):k*Wb
        sum10 = sum10+rt1(n).*s0((n-(k-1)*Wb));
        sum11 = sum11+rt1(n).*s1((n-(k-1)*Wb));
        sum20 = sum20+rt2(n).*s0((n-(k-1)*Wb));
        sum21 = sum21+rt2(n).*s1((n-(k-1)*Wb));
    end
    r10(k)=sum10;
    r11(k)=sum11;
    r20(k)=sum20;
    r21(k)=sum21;
end
% step k is applied also for snr_db = 0
%% step j
bn = 1:1:N;
figure(5)
scatter(bn,r10);
xlabel("b[n]");
hold on
scatter(bn,r11);
xlabel("b[n]");
ylabel("r0[b] and r1[b]");
legend("r0[b] for 15 dB SNR","r1[b] for 15 dB SNR");
title("Correlator Outputs r0[k] and r1[k] with respect to each b[n] for 15 dB SNR");

figure(6)
scatter(bn,r20);
xlabel("b[n]");
hold on
scatter(bn,r21);
xlabel("b[n]");
ylabel("r0[b] and r1[b]");
legend("r0[b] for 0 dB SNR","r1[b] for 0 dB SNR");
title("Correlator Outputs r0[k] and r1[k] with respect to each b[n] for 0 dB SNR");