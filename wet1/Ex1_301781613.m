% Students:
%  Itay Levi - 203192216 
%  Hagai Kariti - 301781613

%% Clear all
clear all;
close all;
clc
%% Q1
% define the signal
x_w = @(w,t) sin(w.*t).*(t >= 0).*exp(-2*t);
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
%% Q2
% define the signal x(t)
w1 = 2*pi*250;
w2 = 2*pi*315;
w3 = 2*pi*375;
w4 = 2*pi*500;
x = @(t) x_w(w1, t).*(t >= 0 & t < 1) ...
    + x_w(w2, t-1).*(t >= 1 & t < 2) ...
    + x_w(w3, t-2).*(t >= 2 & t < 3) ...
    + x_w(w4, t-3).*(t >= 3 & t < 4);
% calculate the signal in the required time interval
t = linspace(0, 4, 64000);
signal = x(t);
% play it
sample_freq = 64000/4;
soundsc(signal, sample_freq);
%% Q3
% calculate time interval
start_index = 1 + sample_freq * 2.99;
end_index = 1 + sample_freq * 3.01;
interval = t(start_index:end_index);
% plot the signal twice (title and legend will be set by later questions)
figure(2);
subplot(2, 1, 1);
plot(interval, signal(start_index:end_index));
xlabel('Time [s]');
ylabel('Amplitude');
xlim([interval(1) interval(end)]);
hold all;

subplot(2, 1, 2);
plot(interval, signal(start_index:end_index));
xlabel('Time [s]');
ylabel('Amplitude');
xlim([interval(1) interval(end)]);
hold all;
%% Q4
% downsample signal
downsampled_frequency = 2000;
downsample_factor = sample_freq/downsampled_frequency;
downsampled_signal = downsample(signal, downsample_factor);
% calculate time interval
start_index = 1 + 2.99 * downsampled_frequency;
end_index = 1 + 3.01 * downsampled_frequency;
downsampled_interval = downsample(interval, downsample_factor);
% plot
subplot(2, 1, 1);
stem(downsampled_interval, downsampled_signal(start_index:end_index));
title("Continunous signal and 2Khz-sampled signal in times [2.99,3.01]s");
legend('Continuous [16KHz]', 'Sampled (2KHz)');
