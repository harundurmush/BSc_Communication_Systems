clc;
clear all;
close all;
t = linspace(0,2*pi) ;
x = sin(pi/4*t) ;
plot(t,x,'r')
hold on
x = sin(pi/4*t+pi/2) ;  % phase shift by pi/2, which is same as cos 
plot(t,x,'b')