close all
clear
clc
EEG_Signal=load("NewEEGSignal.mat");
EEG_Signal=EEG_Signal.NewEEGSignal;
%% Part 1
autocorr_EEG = xcorr(EEG_Signal);
[n, edges] = hist(abs(autocorr_EEG), 100);
PDF = n / sum(n); %normalization
% Plot the autocorrelation
figure;
subplot(2, 1, 1);
plot(autocorr_EEG);
title('Autocorrelation of EEG Signal');
xlabel('Lag');

% Plot the PDF
subplot(2, 1, 2);
bar(edges, PDF, 'hist');
title('PDF of Autocorrelation');
xlabel('Autocorrelation Value');
ylabel('Probability Density');


%% Part 2
% periodgram
[Pxx_periodogram, F_periodogram] = periodogram(EEG_Signal);

% pwelch parameters
window_size = 256; 
overlap = 128;  

% Compute the PSD using pwelch
[Pxx_pwelch, F_pwelch] = pwelch(EEG_Signal, window_size, overlap);

% Calculate the PDF from the PSD
PDF_periodogram = Pxx_periodogram / sum(Pxx_periodogram);
PDF_pwelch = Pxx_pwelch / sum(Pxx_pwelch);

% Plot the PDFs
figure;
subplot(2, 1, 1);
bar(F_periodogram, PDF_periodogram, 'hist');
title('Periodogram');
xlabel('Frequency (Hz)');
ylabel('Probability Density');

subplot(2, 1, 2);
bar(F_pwelch, PDF_pwelch, 'hist');
title('pwelch');
xlabel('Frequency (Hz)');
ylabel('Probability Density');
