close all
clear 
clc
%% Q2 Part 1

% Load EEG signal
EEG_Signal = load('NewEEGSignal.mat');
EEG_Signal = EEG_Signal.NewEEGSignal;
fs = 256;
EEG_Signal_length=length(EEG_Signal);

% Time domain plot
figure;
t = (0:(EEG_Signal_length-1)) / fs;
subplot(3,1,1);
plot(t, EEG_Signal);
title('EEG Signal in the Time Domain');
xlim([0,2])
xlabel('Time (s)');
ylabel('Amplitude');

% Frequency domain plot
f = linspace(0, fs, EEG_Signal_length); 
f_cutoff = 60; 
f_plot = f(f <= f_cutoff);
EEG_Signal_fft = abs(fft(EEG_Signal)) / EEG_Signal_length;
EEG_Signal_fft_plot = EEG_Signal_fft(1:length(f_plot));
subplot(3,1,2);
plot(f_plot, EEG_Signal_fft_plot);
title('EEG Signal in the Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 f_cutoff]); 

% Short-Time Fourier Transform (STFT)
window_length = 128; %Window length
overlap = 64; % Overlap value
[S, F, T, P] = spectrogram(EEG_Signal, hamming(window_length), overlap, window_length, fs);
subplot(3,1,3);
imagesc(T, F, 10*log10(abs(P)));
axis xy;
colormap('jet');
title('STFT of EEG Signal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;

% frequency axis 
ylim([0 f_cutoff]);


%% another way is Power Spectral Density (PSD)
NFFT = 2^nextpow2(EEG_Signal_length);  
f = (0:NFFT/2-1) * (fs / NFFT);  % Frequency axis
EEG_Signal_fft = fft(EEG_Signal, NFFT) / EEG_Signal_length;  
PSD = 2 * abs(EEG_Signal_fft(1:NFFT/2)).^2;  

% Plot the PSD
figure;
plot(f, 10*log10(PSD));
title('Power Spectral Density (PSD) of EEG Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');

% Set a frequency range 
f_cutoff = 60;  
xlim([0 f_cutoff]);

%% Q2 Part2
% frequencis for downsampling
fs_original = 256;  
fs_target = 128;   

% Design a low-pass filter
cutoff_freq = fs_target / 2;  
order = 10;  
lowpass_filter = designfilt('lowpassfir', 'FilterOrder', order, 'CutoffFrequency', cutoff_freq, 'SampleRate', fs_original);

% Filtering
EEG_Signal_filtered = filtfilt(lowpass_filter, EEG_Signal);

% resample factor
downsampling_factor = fs_original / fs_target;

% downsampling
EEG_Signal_ds = downsample(EEG_Signal_filtered, downsampling_factor);

% Time axis
t_downsampled = (0:(length(EEG_Signal_ds) - 1)) / fs_target;

% Time domain plot 
figure;
subplot(3, 1, 1);
plot(t_downsampled, EEG_Signal_ds);
title('Downsampled EEG Signal at 128 Hz (Time Domain)');
xlabel('Time (s)');
ylabel('Amplitude');

% Frequency domain plot 
f_downsampled = linspace(0, fs_target, length(EEG_Signal_ds));
EEG_Signal_ds_fft = abs(fft(EEG_Signal_ds)) / length(EEG_Signal_ds);
subplot(3, 1, 2);
plot(f_downsampled, 10 * log10(EEG_Signal_ds_fft));
title('Downsampled EEG Signal in the Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB/Hz)');

% Short-Time Fourier Transform (STFT) 
window_length = 128;
overlap = 64;
[S_ds, F_ds, T_ds, P_ds] = spectrogram(EEG_Signal_ds, hamming(window_length), overlap, window_length, fs_target);
subplot(3, 1, 3);
imagesc(T_ds, F_ds, 10 * log10(abs(P_ds)));
axis xy;
colormap('jet');
title('STFT of Downsampled EEG Signal');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
colorbar;


%% Q2 Part3
N = 256; % Number of samples

% DFT for EEG_Signal_ds
fs_target = 128; 

% DFT parameters
window_lengths = [N/2, N/4, N/8];
figure;
for i = 1:3
    hamm_win = hamming(window_lengths(i));
    window_length = window_lengths(i);
    
    % Apply the Hamming window to EEG_Signal_ds
    windowed_signal = EEG_Signal_ds(1:window_length)' .* hamm_win;
    
    % Calculate the DFT
    X = fft(windowed_signal, window_length);
    F = fs_target * (0:(window_length - 1)) / window_length;
    
    subplot(3, 1, i);
    plot(F, 10 * log10(abs(X)));
    title(['DFT of Windowed EEG Signal (Hamming Window, ', num2str(N / (2^i)), ' samples)']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    xlim([0, fs_target / 2]);
end


%% Q2 Part4
N = 256; % Number of samples

% DFT parameters
window_lengths = [N/2, N/4, N/8];
figure;
for i = 1:3
    hamm_win = hamming(window_lengths(i));
    window_length = window_lengths(i);
    
    % Apply the Hamming window 
    windowed_signal = EEG_Signal_ds(1:window_length)' .* hamm_win;
    
    % Zeropadding the windowed signal 
    zero_arr=zeros([1,N]);
    zero_arr(1:window_length)=windowed_signal;
     padded_signal = zero_arr;
    % Calculate the N-point DFT
    X = fft(padded_signal, N);
    F = fs_target * (0:(N - 1)) / N;
    
    subplot(3, 1, i);
    plot(F, 10 * log10(abs(X)));
    title(['N-Point DFT of Windowed and Zero-Padded EEG Signal (Hamming Window, ', num2str(N / (2^i)), ' samples)']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    xlim([0, fs_target / 2]);
end
