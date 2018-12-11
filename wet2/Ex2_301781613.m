% Students:
%  Itay Levi - 203192216 
%  Hagai Kariti - 301781613
%% Clear all
clear all; close all; clc
%% Q1
% run fft on input signal
load('signal.mat');
fft_x = fft(x);
sample_rate = 2000;
% plot input signal (positive freq only)
f_axis = linspace(0, sample_rate/2, length(x)/2);
positive_fft_x = fft_x(1:length(f_axis));
figure(1);
subplot(2, 1, 1);
plot(f_axis, mag2db(abs(positive_fft_x)));
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
figure(4);
subplot(4, 1, 1);
plot(f_axis, mag2db(abs(positive_fft_x)));
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
%% Q3
function fft_matrix = fft_windows(X, N)
    length_remainder = mod(length(X), N);
    padded_X = padarray(X, [0, length_remainder ], 'post');
    windows_mat = reshape(padded_X, length(padded_X)/N, N);
    fft_matrix = fft(windows_mat);
end