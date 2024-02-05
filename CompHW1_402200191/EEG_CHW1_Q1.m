close all
clear
clc
%%  Q1 Part1
t=0:1e-3:2;
x=chirp(t,100,2,500,'quadratic');
plot(x)
%% Q1 Part2
N=128; % length
n=0:N-1; %samples
% Windows
Rect_Win = rectwin(N);
Tri_Win = triang(N);
Gauss_Win = gausswin(N);
Hamm_Win = hamming(N);

wvtool(Rect_Win);
wvtool(Tri_Win);
wvtool(Gauss_Win);
wvtool(Hamm_Win);
 

%% Q1 Part 3
[S_R, f_R, t_R] = spectrogram(x, Rect_Win, 0, 128, 1000);

% Plot the spectrograms

figure;
subplot(2,2,1)
imagesc(t_R, f_R, 10*log10(abs(S_R)));
axis xy;
colormap('jet');
title('Spectrogram with Rectangular Window');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;
[S_Tri, f_Tri, t_Tri] = spectrogram(x, Tri_Win, 0, 128, 1000);

subplot(2,2,2)
imagesc(t_Tri, f_Tri, 10*log10(abs(S_Tri)));
axis xy;
colormap('jet');
title('Spectrogram with Triangular Window');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

[S_Gauss, f_Gauss, t_Gauss] = spectrogram(x, Gauss_Win, 0, 128, 1000);

subplot(2,2,3)
imagesc(t_Gauss, f_Gauss, 10*log10(abs(S_Gauss)));
axis xy;
colormap('jet');
title('Spectrogram with Gaussian Window');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

[S_Hamm, f_Hamm, t_Hamm] = spectrogram(x, Hamm_Win, 0, 128, 1000);

subplot(2,2,4)
imagesc(t_Hamm, f_Hamm, 10*log10(abs(S_Hamm)));
axis xy;
colormap('jet');
title('Spectrogram with Hamming Window');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;


%% Q1 Part4
% I chose the Hamming Widow 
[S_Hamm0, f_Hamm0, t_Hamm0] = spectrogram(x, Hamm_Win, 0, 128, 1000);
[S_Hamm32, f_Hamm32, t_Hamm32] = spectrogram(x, Hamm_Win, 32, 128, 1000);
[S_Hamm64, f_Hamm64, t_Hamm64] = spectrogram(x, Hamm_Win, 64, 128, 1000);
[S_Hamm127, f_Hamm127, t_Hamm127] = spectrogram(x, Hamm_Win, 127, 128, 1000);
figure;
% Plot the spectrograms
% No overlap
subplot(2,2,1);
imagesc(t_Hamm0, f_Hamm0, 10*log10(abs(S_Hamm0)));
axis xy;
colormap('jet');
title('Hamming Window (No Overlap)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Overlap 32
subplot(2,2,2);
imagesc(t_Hamm32, f_Hamm32, 10*log10(abs(S_Hamm32)));
axis xy;
colormap('jet');
title('Hamming Window (Overlap: 32)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Overlap 64
subplot(2,2,3);
imagesc(t_Hamm64, f_Hamm64, 10*log10(abs(S_Hamm64)));
axis xy;
colormap('jet');
title('Hamming Window (Overlap: 64)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Overlap 127
subplot(2,2,4);
imagesc(t_Hamm127, f_Hamm127, 10*log10(abs(S_Hamm127)));
axis xy;
colormap('jet');
title('Hamming Window (Overlap: 127)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

sgtitle('Spectrogram Comparison for Hamming Window with Different Overlaps');


%% Q1 Part 5
% Windows
Hamm_Win32 = hamming(32);
Hamm_Win64 = hamming(64);
Hamm_Win128 = hamming(128);
Hamm_Win512 = hamming(512);
% Spectrograms
[S_Hamm32L, f_Hamm32L, t_Hamm32L] = spectrogram(x, Hamm_Win32, 31, 32, 1000);
[S_Hamm64L, f_Hamm64L, t_Hamm64L] = spectrogram(x, Hamm_Win64, 63, 64, 1000);
[S_Hamm128L, f_Hamm128L, t_Hamm128L] = spectrogram(x, Hamm_Win128, 127, 128, 1000);
[S_Hamm512L, f_Hamm512L, t_Hamm512L] = spectrogram(x, Hamm_Win512, 511, 512, 1000);

% Plot the spectrograms
figure;

% Window Length 32
subplot(2,2,1);
imagesc(t_Hamm32L, f_Hamm32L, 10*log10(abs(S_Hamm32L)));
axis xy;
colormap('jet');
title('Hamming Window (L=32, Overlap=31)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Window Length 64
subplot(2,2,2);
imagesc(t_Hamm64L, f_Hamm64L, 10*log10(abs(S_Hamm64L)));
axis xy;
colormap('jet');
title('Hamming Window (L=64, Overlap=63)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Window Length 128
subplot(2,2,3);
imagesc(t_Hamm128L, f_Hamm128L, 10*log10(abs(S_Hamm128L)));
axis xy;
colormap('jet');
title('Hamming Window (L=128, Overlap=127)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% Window Length 512
subplot(2,2,4);
imagesc(t_Hamm512L, f_Hamm512L, 10*log10(abs(S_Hamm512L)));
axis xy;
colormap('jet');
title('Hamming Window (L=512, Overlap=511)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

sgtitle('Spectrogram Comparison for Hamming Window with Different Window Lengths');


%% Q1 Part6
Hamm_Win256=hamming(256);
Hamm_Win1024=hamming(1024);
[S_HammL, f_HammL, t_HammL] = spectrogram(x, Hamm_Win, 64, 128, 1000);
[S_Hamm2L, f_Hamm2L, t_Hamm2L] = spectrogram(x, Hamm_Win256, 128, 256, 1000);
[S_Hamm4L, f_Hamm4L, t_Hamm4L] = spectrogram(x, Hamm_Win512, 256, 512, 1000);
[S_Hamm8L, f_Hamm8L, t_Hamm8L] = spectrogram(x, Hamm_Win1024, 512, 1024, 1000);

% Plot the spectrograms
figure;

% 128
subplot(2,2,1);
imagesc(t_HammL, f_HammL, 10*log10(abs(S_HammL)));
axis xy;
colormap('jet');
title('Hamming Window (Overlap: 64, FFT Length: 128)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% 256
subplot(2,2,2);
imagesc(t_Hamm2L, f_Hamm2L, 10*log10(abs(S_Hamm2L)));
axis xy;
colormap('jet');
title('Hamming Window (Overlap: 128, FFT Length: 256)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% 512
subplot(2,2,3);
imagesc(t_Hamm4L, f_Hamm4L, 10*log10(abs(S_Hamm4L)));
axis xy;
colormap('jet');
title('Hamming Window (Overlap: 256, FFT Length: 512)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% 1024
subplot(2,2,4);
imagesc(t_Hamm8L, f_Hamm8L, 10*log10(abs(S_Hamm8L)));
axis xy;
colormap('jet');
title('Hamming Window (Overlap: 512, FFT Length: 1024)');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

sgtitle('Spectrogram Comparison for Hamming Window with Different Overlaps and FFT Lengths');
