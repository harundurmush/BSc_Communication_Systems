clc;
clear all;
%% 2.1 step a
fs=10;
T=25;
t=0:1/fs:T;
x1=linspace(1,(fs*T)+1);
for k=1:(fs*T)+1
    if k<=fs*5
        x1(k)=2*t(k)-5;
    elseif k>fs*5 && k<fs*15
        x1(k)=5*cos(2*pi*t(k));
    else
        x1(k)=100*exp(-t(k)/5);
    end
end
%% 2.1 step b
figure;
plot(t,x1,"blue");
title("x1(t) signal in time domain");
ylabel("x1(t)");
xlabel("time (s)");
%% 2.1 step c
x2=linspace(1,(fs*T)+1);
for k=1:(fs*T)+1
    x2(k)=cos(4*pi*t(k));
end

figure;
plot(t,x2,"black");
title("**EXTRA** x2(t) signal in time domain");
ylabel("x2(t)");
xlabel("time (s)");

y1=x1.*x2;

figure;
plot(t,y1,"red");
title("y1(t) signal in time domain");
ylabel("y1(t)");
xlabel("time (s)");
%% length reconfiguration of x1 and x2
N=length(x1);
z_matrix=zeros(1,(N-1)/2);
extendedx1= [ z_matrix x1 z_matrix];
extendedx2= [ z_matrix x2 z_matrix];
%% 2.2 step a
N=length(extendedx1); % N is changed %
x1f=fftshift(abs(fft(extendedx1,N)))/N;
y1f_1=fftshift(abs(fft(y1,N)))/N;

figure;
subplot(2,1,1)
plot(x1f,"blue");
title("Fourier Transform of the x1(t) and y1(t) signals");
ylabel("|x1(f)|");
subplot(2,1,2)
plot(y1f_1,"red");
ylabel("|y1(f)|");
xlabel("frequency (Hz)");
%% 2.2 step b
x2f=fftshift(abs(fft(extendedx2,N)))/N;
y2f=x1f.*x2f;

figure;
plot(y2f,"red");
title("Fourier Transform of output signal y2(t)");
ylabel("|y2(f)|");
xlabel("frequency (Hz)");
%% 2.2 step c
y2t=ifft(y2f,N);

figure;
plot(y2t,"blue");
title("Inverse Fourier Transform of output signal y2(f)");
ylabel("y2(t)");
xlabel("time (s)");
%% 2.2 step d

y2t_2=conv(x1,x2);

figure;
plot(y2t_2,"blue");
title("Convulution of input signal and impulse response");
ylabel("y2(t)");
xlabel("time (s)");