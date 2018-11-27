% Students:
%  Itay Levi - 203192216 
%  Hagai Kariti - 301781613

%% Clear all
clear all;
close all;
clc
%% Q1
% define fourier transform for X
w0 = 2*pi*250;
X_fourier = @(w) w0./(w0^2 + (j.*w + 2).^2);
% calculate it for given interval
f_space = linspace(-300, 300, 10000);
X_f_abs = abs(X_fourier(2*pi*f_space));
% plot it
figure(1);
subplot(3, 1, 1);
plot(f_space, X_f_abs);
title("|F\{X_{\omega_0}\}| in frequencies -300 to 300, \omega_0 = 2\pi\cdot 250");
xlabel("Frequency [Hz]");
ylabel("Magnitude");