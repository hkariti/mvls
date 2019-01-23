% Students:
%  Itay Levi - 203192216 
%  Hagai Kariti - 301781613
%% Clear all
clear all; close all; clc
%% Q3
[channel, Fs] = audioread ('channel.wav');
N_4 = 1000;
sample_1_1k=1:N_4;
sample_50k_51k=50001:51000;
sample_100k_101k=100001:101000;
f_4=linspace(-Fs/2+Fs/(2*N_4),Fs/2,N_4);
sample_50k_51k_dft=fft(channel(sample_50k_51k)');
sample_1_1k_dft=fft(channel(sample_1_1k)');
sample_100k_101k_dft=fft(channel(sample_100k_101k)');
figure(1)
subplot(3,1,1)
plot(f_4,abs(fftshift(sample_1_1k_dft)));
title("1-1000 samples of channel");
xlabel("f[Hz]");
ylabel("|DFT(channel)|");
subplot(3,1,2)
plot(f_4,abs(fftshift(sample_50k_51k_dft)));
title("50k-51k samples of channel");
xlabel("f[Hz]");
ylabel("|DFT(channel)|");
subplot(3,1,3)
plot(f_4,abs(fftshift(sample_100k_101k_dft)));
title("100k-101k samples of channel");
xlabel("f[Hz]");
ylabel("|DFT(channel)|");
%% temporary
t = 0:(1/Fs):(length(channel)-1)/Fs;
x_AM = [ t',channel];
fc1 = 1e4;
fc2 = 2e4;
fc3 = 3e4;
%BW = 8e3;
%% Q5
x = [((0:(length(channel)-1))/Fs)',channel];
coeff = load('coeff.mat');
coeff = coeff.coeff;
figure(2);
set(gcf,'color','w');
stem(coeff); 
set(xlabel('$index$'),'interpreter','Latex');
set(title('$Coefficients$'),'interpreter','Latex'); 
xlim([0 length(coeff)+1]);
%% Q6
figure(3);
freqz(coeff,1,[],Fs);
set(gcf,'color','w');
set(title('Amplitude and phase response'));
%% Q7
sim('simulink_sim');
x_AM_filter = (x_AM_filt.Data);
dft_x_AM_filter = fftshift(fft(x_AM_filter(100001:101000)));
figure(4)
subplot(2,1,1)
plot(f_4,abs(dft_x_AM_filter));
title('DFT(x_A_M filter)');
xlabel('f[Hz]');
ylabel('|DFT(x_A_M filter)|');
%% Q8
x_AM_modulation = (x_AM_mod.Data);
dft_x_AM_modulation = fftshift(fft(x_AM_modulation(100001:101000)));
subplot(2,1,2)
plot(f_4,abs(dft_x_AM_modulation));
title('DFT(x_A_M modulation)');
xlabel('f[Hz]');
ylabel('|DFT(x_A_M modulation)|');
%% Q11
N1=60;
n1=0:N1;
x1=sinc(n1-N1/2);
x2=0.625*sinc(0.625*(n1-N1/2));
h1=(x1-x2).*blackman(N1+1)';
figure(5);
subplot(2,1,2);
freqz(h1,1,[],Fs);
set(gcf,'color','w');
set(title('IIR - phase and amplitude response')); 
%% Q12
N2=130;
n2=0:N2;
y1=sinc(n2-N2/2);
y2=0.625*sinc(0.625*(n2-N2/2));
h2=(y1-y2).*blackman(N2+1)';
figure(6); 
freqz(h2,1,[],Fs);
set(gcf,'color','w');
set(title('FIR - phase and amplitude response'));