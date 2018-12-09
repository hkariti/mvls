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
sample_frequency_continuous = 64000/4;
soundsc(signal, sample_frequency_continuous);
%% Q3
% calculate time interval
start_index = 1 + sample_frequency_continuous * 2.99;
end_index = 1 + sample_frequency_continuous * 3.01;
interval = t(start_index:end_index);
% plot the signal twice (title and legend will be set by later questions)
figure(2);
subplot(2, 1, 1);
plot(interval, signal(start_index:end_index));
xlabel('Time [s]');
ylabel('Amp');
xlim([interval(1) interval(end)]);
hold all;

subplot(2, 1, 2);
plot(interval, signal(start_index:end_index));
xlabel('Time[sec]');
ylabel('Amplitude');
xlim([interval(1) interval(end)]);
hold all;
%% Q4
% downsample signal
sample_frequency = 2000;
downsample_factor = sample_frequency_continuous/sample_frequency;
sampled_signal = downsample(signal, downsample_factor);
% calculate time interval
start_index = 1 + 2.99 * sample_frequency;
end_index = 1 + 3.01 * sample_frequency;
sampled_interval = downsample(interval, downsample_factor);
% plot
subplot(2, 1, 1);
stem(sampled_interval, sampled_signal(start_index:end_index));
title(['Original signal v.s. Reconstructions, Fs = ', num2str(sample_frequency), '[Hz]']);
legend('Continuous [16KHz]', 'Sampled (2KHz)');
%% Q5
% Calculate fft on sampled signal
signal_fft = fftshift(fft(sampled_signal));
% plot it
num_of_samples = length(sampled_signal);
frequencies = linspace(-sample_frequency/2, sample_frequency/2, num_of_samples);
figure(1);
subplot(3, 1, 2);
plot(frequencies, abs(signal_fft));
title('FFT of 2KHz sampled signal');
xlabel('Frequencies [Hz]');
ylabel('Magnitude');
%% Q7
x_up1 = upsample(sampled_signal,downsample_factor);
% t=t2;

m_US1 = find((t>=2.99) & (t<=3.01));
dt = t(m_US1);

Int_sinc=sinc((-256:255)/downsample_factor);
x_ID1_conv=conv(x_up1,Int_sinc);
x_ID1=x_ID1_conv(256:(length(x_ID1_conv)-256));

figure(2);
subplot(2,1,1);
plot(dt,x_ID1(m_US1));
hold all;
%% Q8
Int_ZOH1=linspace(1,1,downsample_factor);
x_ZOH1_conv=conv(x_up1,Int_ZOH1);
x_ZOH1=x_ZOH1_conv(1:(length(x_ZOH1_conv)-downsample_factor+1));

figure(2);
subplot(2,1,1);
plot(dt,x_ZOH1(m_US1));
hold all;
%% Q9
Int_FOH1=conv(Int_ZOH1,Int_ZOH1)/downsample_factor;
x_FOH1_conv=conv(x_up1,Int_FOH1);
x_FOH1=x_FOH1_conv(downsample_factor:(length(x_FOH1_conv)-downsample_factor+1));

figure(2);
subplot(2,1,1);
plot(dt,x_FOH1(m_US1));
hold all;

legend('original','sampled','ideal','ZOH','FOH');
%% Q11
Fs2=800;
N_DS2=sample_frequency_continuous/Fs2;
x_DS2=downsample(signal,N_DS2);
t_DS2=downsample(t,N_DS2);

m11=find((t_DS2>=2.99) & (t_DS2<=3.01));
dx_DS2=x_DS2(m11);
dt_DS2=t_DS2(m11);

N2=4*Fs2;
x11=fft(x_DS2);
w2_FFT=linspace(-Fs2/2,Fs2/2,N2);
x2_FFT = fftshift( x11 );

figure(1);
subplot(3,1,3);
plot(w2_FFT,abs(x2_FFT));
title(['F_s = ',num2str(Fs2), 'Hz']);
xlabel('f[Hz]');
ylabel('|DFT(X)|');
hold all;


figure(2);
subplot(2,1,2);
stem(dt_DS2,dx_DS2)
title(['F_s = ',num2str(Fs2), 'Hz']);
xlabel('Time[sec]');
ylabel('Amp');
hold all;

x_up2 = upsample(x_DS2,N_DS2);
t_US2 = t;

m_US2=find((t_US2>=2.99) & (t_US2<=3.01));
dt_US2=t_US2(m_US2);

Int_sinc=sinc((-256:255)/N_DS2);
x_ID2_conv=conv(x_up2,Int_sinc);
x_ID2=x_ID2_conv(256:(length(x_ID2_conv)-256));

figure(2);
subplot(2,1,2);
plot(dt_US2,x_ID2(m_US2));
hold all;

Int_ZOH2=linspace(1,1,N_DS2);
x_ZOH2_conv=conv(x_up2,Int_ZOH2);
x_ZOH2=x_ZOH2_conv(1:(length(x_ZOH2_conv)-N_DS2+1));

figure(2);
subplot(2,1,2);
plot(dt_US2,x_ZOH2(m_US2));
hold all;

Int_FOH2=conv(Int_ZOH2,Int_ZOH2)/N_DS2;
x_FOH2_conv=conv(x_up2,Int_FOH2);
x_FOH2=x_FOH2_conv(N_DS2:(length(x_FOH2_conv)-N_DS2+1));

figure(2);
subplot(2,1,2);
plot(dt_US2,x_FOH2(m_US2));
hold all;

legend('original','sampled','ideal','ZOH','FOH');