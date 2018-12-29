% Students:
%  Itay Levi - 203192216 
%  Hagai Kariti - 301781613
%% Clear all
clear all; close all; clc
%% Q1
% run fft on input signal
load('signal.mat');
fft_x = fft(x);
x_len = length(fft_x);
sample_rate = 2000;
% plot input signal (positive freq only)
f_axis = linspace(0, sample_rate/2, length(x)/2);
positive_fft_x = fft_x(1:length(f_axis));
figure(1);
subplot(2, 1, 1);
plot(f_axis, mag2db(abs(positive_fft_x)));
title ('DFT of x[n]');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
figure(4);
subplot(4, 1, 1);
plot(f_axis, mag2db(abs(positive_fft_x)));
title ('DFT of x[n]');
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
%% q4
N_4=256;
F_4=2000;
T_4=1/F_4;
M_4=ceil(x_len/N_4);
t_total_4=length(x)*T_4;
t_4=linspace(0,t_total_4,M_4);
FFT_x_4=fft_windows(x,N_4);
f_4=linspace(0,F_4/2,N_4/2);
figure(1)
subplot(2,1,2);
mesh(f_4,t_4,20*log10(abs(FFT_x_4(:,1:N_4/2))));
title('3D spectrogram');
xlabel('Frequency[Hz]');
ylabel('Time[sec]');
zlabel('DFT-X [dB]');
view(15,75);

%% q5
t_find_5=0.55;
interval_t_5=t_total_4/M_4;
n_5=ceil(t_find_5/interval_t_5);
DFT_find_5=abs(FFT_x_4(n_5,:));
figure(2)
subplot(2,3,1);
plot(f_4,20*log10(DFT_find_5(1:N_4/2)));
title(['DFT-X between ',num2str((n_5-1)*0.128),'[sec] to ',num2str((n_5)*0.128),'[sec]']);
xlabel('Frequency[Hz]');
ylabel('DFT-X [dB]');
hold on;
%% q6
x_5=x(((n_5-1)*N_4+1):(n_5*N_4));
window_length_6=128;
x_6=x_5((N_4/2-window_length_6/2+1):(N_4/2+window_length_6/2));
DFT_6=abs(fft(x_6));
f_6=linspace(0,F_4/2,window_length_6/2);
subplot(2,3,2)
plot(f_6,20*log10(DFT_6(1:window_length_6/2)));
title('DFT-X(128bit unit window)');
xlabel('Frequency[Hz]');
ylabel('DFT-X [dB]');
%% q7
window_7=zeros(1,N_4);
window_7((N_4/2-window_length_6/2+1):(N_4/2+window_length_6/2))=ones(1,window_length_6);
x_7=x_5.*window_7;
DFT_7=abs(fft(x_7));
subplot(2,3,3);
plot(f_4,20*log10(DFT_7(1:N_4/2)));
title('DFT-X(128bit unit window)z.padding');
xlabel('Frequency[Hz]');
ylabel('DFT-X [dB]');

%% q8.5
subplot(2,3,4);
plot(f_4,20*log10(DFT_find_5(1:N_4/2)));
title('DFT-X, fifth interval');
xlabel('Frequency[Hz]');
ylabel('DFT-X [dB]');


%% q8.6+q8.7
N_8=window_length_6;
n_8=0:(N_8-1);
blackman_8=0.42-0.5*cos(2*pi*n_8/(N_8-1))+0.08*cos(4*pi*n_8/(N_8-1));
window_8_7=zeros(1,N_4);
window_8_7((N_4/2-N_8/2+1):(N_4/2+N_8/2))=blackman_8;
windowed_x_8_6=x_6.*blackman_8;
windowed_x_8_7 =x_5.*window_8_7;
DFT_x_8_6=fft(windowed_x_8_6);
DFT_x_8_7=fft(windowed_x_8_7);
subplot(2,3,5)
plot(f_6,20*log10(abs(DFT_x_8_6(1:N_8/2))));
title('DFT-X - Blackman no z.padding');
xlabel('Frequency[Hz]');
ylabel('DFT-X [dB]');
subplot(2,3,6)
plot(f_4,20*log10(abs(DFT_x_8_7(1:N_4/2))));
title('DFT-X - Blackman z.padding');
xlabel('Frequency[Hz]');
ylabel('DFT-X [dB]');
%% Question 12
Filter_1 = [1,1];
Filter_2 = [1,-1];
Filter_3 = [1,0,0,-1];

figure(3)
hold on
subplot(1,3,1)
zplane(Filter_1)
title ('zeros&poles map for Filter_1')
subplot(1,3,2)
zplane(Filter_2)
title ('zeros&poles map for Filter_2')
subplot(1,3,3)
zplane(Filter_3)
title ('zeros&poles map for Filter_3')

    %% Question 13
    x_Filter_1 = conv(x,Filter_1);
    x_Filter_2 = conv(x,Filter_2);
    x_Filter_3 = conv(x,Filter_3);

    f_Filter_1 = linspace(-sample_rate/2,sample_rate/2,length(x_Filter_1));
    f_Filter_2 = linspace(-sample_rate/2,sample_rate/2,length(x_Filter_2));
    f_Filter_3 = linspace(-sample_rate/2,sample_rate/2,length(x_Filter_3));

    x_dft_Filter_1 = fftshift(fft(x_Filter_1));
    x_dft_Filter_2 = fftshift(fft(x_Filter_2));
    x_dft_Filter_3 = fftshift(fft(x_Filter_3));

    figure(4)
    subplot(4,1,2)
    plot(f_Filter_1, 20*log10(abs(x_dft_Filter_1)))
    xlim([1,sample_rate/2])
    grid on
    title('DFT of x[n] filtered by Filter_1')
    xlabel('Freq[Hz]')
    ylabel('Magnitude [dB]')
    subplot(4,1,3)
    plot(f_Filter_2, 20*log10(abs(x_dft_Filter_2)))
    xlim([1,sample_rate/2])
    grid on
    title('DFT of x[n] filtered by Filter_2')
    xlabel('Freq[Hz]')
    ylabel('Magnitude [dB]');
    subplot(4,1,4)
    plot(f_Filter_3, 20*log10(abs(x_dft_Filter_3)))
    xlim([1,sample_rate/2])
    grid on
    title('DFT of x[n] filtered by Filter_3')
    xlabel('Freq[Hz]')
    ylabel('Magnitude [dB]');