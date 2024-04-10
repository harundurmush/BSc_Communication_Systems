clc;
clear all;
close all;
%% 2.1 DM Modulator
%% 2.1.a
Td = 0.05;
Fs = 2000;
Ts = 1/Fs;
time = 0:Ts:(Td-Ts);
%% 2.1.b
f1 = 60;
f2 = 90;
mt = 0.8*cos(2*pi*f1.*time) + 0.7*sin(2*pi*f2.*time);
%% 2.1.c Report
%% 2.1.d
n = length(mt);
delta2 = 0.1;
delta3 = 0.5;
for i = 1:(n-1) 
    m_der(i) = (mt(i+1) - mt(i))/Ts;
end
delta_min=max(abs(m_der))/Fs;
me1(1)=0;
err1(1)=0;
errq1(1)=0;
for j=2:n
    err1(j)=mt(j)-me1(j-1);
    errq1(j)=delta_min*sign(err1(j));
    me1(j)=me1(j-1)+errq1(j);
end
me2(1)=0;
err2(1)=0;
errq2(1)=0;
for j=2:n
    err2(j)=mt(j)-me2(j-1);
    errq2(j)=delta2*sign(err2(j));
    me2(j)=me2(j-1)+errq2(j);
end
me3(1)=0;
err3(1)=0;
errq3(1)=0;
for j=2:n
    err3(j)=mt(j)-me3(j-1);
    errq3(j)=delta3*sign(err3(j));
    me3(j)=me3(j-1)+errq3(j);
end
figure (4)
subplot(311)
plot(time,mt);
ylabel("original signal m(t)");
hold on
stairs(time,me1);
xlabel("time (s)");
ylabel("modulated signal me1(t)");
legend("m(t)","me1(t) OPTIMUM delta");
title("Staircase Approximation With Different Delta Values");

subplot(312)
plot(time,mt);
ylabel("original signal m(t)");
hold on
stairs(time,me2);
xlabel("time (s)");
ylabel("modulated signal me2(t)");
legend("m(t)","me2(t) delta=0.1");

subplot(313)
plot(time,mt);
ylabel("original signal m(t)");
hold on
stairs(time,me3);
xlabel("time (s)");
ylabel("modulated signal me3(t)");
legend("m(t)","me3(t) delta=0.5");
encoder1(1) = 0;
for k=2:n
    if(sign(err1(k)) == 1)
        encoder1(k) = 1;
    else
        encoder1(k) = 0;
    end
end
encoder2(1) = 0;
for k=2:n
    if(sign(err2(k)) == 1)
        encoder2(k) = 1;
    else
        encoder2(k) = 0;
    end
end
encoder3(1) = 0;
for k=2:n
    if(sign(err3(k)) == 1)
        encoder3(k) = 1;
    else
        encoder3(k) = 0;
    end
end
%% 2.2 DM Demodulator
%%
n1 = length(me1);
n2 = length(me2);
n3 = length(me3);
md1(1) = 0;
for i=2:n1
    md1(i) = me1(i);
    if(encoder1(i) ~= 0)
    md1(i) = md1(i-1) + delta_min;
    elseif(encoder1(i) == 0)
    md1(i) = md1(i-1) - delta_min;
    end
end
md2(1) = 0;
for i=2:n2
    md2(i) = me2(i);
    if(encoder2(i) ~= 0)
    md2(i) = md2(i-1) + delta2;
    elseif(encoder2(i) == 0)
    md2(i) = md2(i-1) - delta2;
    end
end
md3(1) = 0;
for i=2:n3
    md3(i) = me3(i);
    if(encoder3(i) ~= 0)
    md3(i) = md3(i-1) + delta3;
    elseif(encoder3(i) == 0)
    md3(i) = md3(i-1) - delta3;
    end
end
%% low-pass filtering
%%
cf=0.1;
n=100;
b=fir1(n,cf);
out1=conv2(md1,b,'same');
out2=conv2(md2,b,'same');
out3=conv2(md3,b,'same');
%% 3 Results
%%
figure (1)
subplot(211)
plot(time,mt);
xlabel("time (s)");
ylabel("m(t)");
hold on
plot(time,md1);
xlabel("time (s)");
ylabel("mdem(t)");
legend("m(t)","mdem(t)");
title("Comparison between original signal m(t) and demodulated signal mdem(t) with OPTIMUM delta");

subplot(212)
plot(time,mt);
xlabel("time (s)");
ylabel("m(t)");
hold on
plot(time,out1);
xlabel("time (s)");
ylabel("mrec(t)");
legend("m(t)","mrec(t)");
title("Comparison between original signal m(t) and reconstructed signal mrec(t) with OPTIMUM delta");

figure(2)
subplot(211)
plot(time,mt);
xlabel("time (s)");
ylabel("m(t)");
hold on
plot(time,md2);
xlabel("time (s)");
ylabel("mdem(t)");
legend("m(t)","mdem(t)");
title("Comparison between original signal m(t) and demodulated signal mdem(t) with delta = 0.1");

subplot(212)
plot(time,mt);
xlabel("time (s)");
ylabel("m(t)");
hold on
plot(time,out2);
xlabel("time (s)");
ylabel("mrec(t)");
legend("m(t)","mrec(t)");
title("Comparison between original signal m(t) and reconstructed signal mrec(t) with delta = 0.1");

figure(3)
subplot(211)
plot(time,mt);
xlabel("time (s)");
ylabel("m(t)");
hold on
plot(time,md3);
xlabel("time (s)");
ylabel("mdem(t)");
legend("m(t)","mdem(t)");
title("Comparison between original signal m(t) and demodulated signal mdem(t) with delta = 0.5");

subplot(212)
plot(time,mt);
xlabel("time (s)");
ylabel("m(t)");
hold on
plot(time,out3);
xlabel("time (s)");
ylabel("mrec(t)");
legend("m(t)","mrec(t)");
title("Comparison between original signal m(t) and reconstructed signal mrec(t) with delta = 0.5");
