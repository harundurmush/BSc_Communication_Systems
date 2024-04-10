clc;
clear all;
close all;
%% Creating message signal
%% step 1
Fs = 27*10^3;
Rb = 18*10^3;
Ts = 1/Fs;
Tb = 1/Rb;
N = 1000;
M = 4;
W = 6000;
A = 1;
randoms = randi(4,1,N);
k = log2(M);
%% step 2
for j=1:N
    if(randoms(j)==1)
        pamd(j)=-3*A;
    elseif(randoms(j)==2)
        pamd(j)=-1*A;
    elseif(randoms(j)==3)
        pamd(j)=1*A;
    elseif(randoms(j)==4)
        pamd(j)=3*A;
    end
end
%% step 3
Rs = Rb/k;
sps = Fs/Rs;
%% Pulse shaping the data
%% step 4
rolloff_1 = 0;
rolloff_2 = 0.5;
rolloff_3 = 1;
span = 10;

ht0 = rcosdesign(rolloff_1,span,sps,"sqrt");
ht1 = rcosdesign(rolloff_2,span,sps,"sqrt");
ht2 = rcosdesign(rolloff_3,span,sps,"sqrt");

L = length(ht0);

t = linspace(-N/Rs,N/Rs,L);

hf0 = abs(fftshift((fft(ht0,L))/L));
hf1 = abs(fftshift((fft(ht1,L))/L));
hf2 = abs(fftshift((fft(ht2,L))/L));

f=linspace(-Fs/2,Fs/2,L);
%% step 5
figure;

subplot(231)
plot(t,ht0); 
title("SR-RCF signal in time domain"); 
xlabel('time (s)'); 
ylabel('amplitude');
legend("roll-off = 0");

subplot(232)
plot(t,ht1); 
title("SR-RCF signal in time domain"); 
xlabel('time (s)'); 
ylabel('amplitude');
legend("roll-off = 0.5");

subplot(233)
plot(t,ht2); 
title("SR-RCF signal in time domain"); 
xlabel("time (s)");
ylabel('amplitude');
legend("roll-off = 1");

subplot(234)
plot(f,hf0); 
title('SR-RCF signal in frequency domain'); 
xlabel('f (Hz)'); 
ylabel('magnitude');
legend("roll-off = 0");

subplot(235)
plot(f,hf1); 
title('SR-RCF signal in frequency domain'); 
xlabel('f (Hz)'); 
ylabel('magnitude');
legend("roll-off = 0.5");

subplot(236)
plot(f,hf2); 
title('SR-RCF signal in frequency domain'); 
xlabel('f (Hz)'); 
ylabel('magnitude');
legend("roll-off = 1");
%% The transmitted/received signal
%%
cf = 1; % no need to filter channel
%% step 6
t1=upfirdn(pamd,ht0,sps,1);
t2=upfirdn(pamd,ht1,sps,1);
t3=upfirdn(pamd,ht2,sps,1);
%% step 7
r1 = awgn(t1,20,'measured');
r2 = awgn(t2,20,'measured');
r3 = awgn(t3,20,'measured');
%% step 8
yout1 = upfirdn(r1,ht0,1,sps);
yout2 = upfirdn(r2,ht1,1,sps);
yout3 = upfirdn(r3,ht2,1,sps);
%% step 9
yout1 = yout1(10+1:end-10);
yout2 = yout2(10+1:end-10);
yout3 = yout3(10+1:end-10);
%% Eye Diagrams
%% step 10
eyediagram(yout1,sps,Ts,0);
xlabel('time (s)'); 
ylabel('amplitude');
title("Eye Diagram for roll-off factor 0");
eyediagram(yout2,sps,Ts,0);
xlabel('time (s)'); 
ylabel('amplitude');
title("Eye Diagram for roll-off factor 0.5");
eyediagram(yout3,sps,Ts,0);
xlabel('time (s)'); 
ylabel('amplitude');
title("Eye Diagram for roll-off factor 1");
%% step 11
rolloff = [0.25 0.33 0.5 0.67 0.75 1];
Rs_v(1) = (2*W)/(1+rolloff(1));
Rs_v(2) = (2*W)/(1+rolloff(2));
Rs_v(3) = (2*W)/(1+rolloff(3));
Rs_v(4) = (2*W)/(1+rolloff(4));
Rs_v(5) = (2*W)/(1+rolloff(5));
Rs_v(6) = (2*W)/(1+rolloff(6));

figure;
plot(rolloff,Rs_v);
grid on
xlabel("roll-off factor values (beta)");
ylabel("symbol rate (Rs)");
title("Symbol rate values with respect to the corresponding roll-off factors");



