clear all
close all
clc
%% Part 1
EEG_Data=load("Ex1.mat").EEG_Sig;
fs=200;
% Define time vector
t = (0:(size(EEG_Data, 2)-1)) / fs;

% Plot EEG data for each channel
figure;

for channel = 1:size(EEG_Data, 1)
    subplot(size(EEG_Data, 1), 1, channel);
    plot(t, EEG_Data(channel, :));
    title(['Channel ' num2str(channel)]);
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    xlim([t(1) t(end)]); 
end
%% part 2
figure;
scatter3(EEG_Data(1, :), EEG_Data(2, :), EEG_Data(3, :), 20, t, 'filled');
colormap('winter');
colorbar;
title('3D Scatter Plot of EEG Data');
xlabel('EEG Channel 1');
ylabel('EEG Channel 2');
zlabel('EEG Channel 3');
grid on;
%% Part 3
EEG_Data_C=EEG_Data-mean(EEG_Data(:));
cov_matrix = cov(EEG_Data_C');
[eigenvectors, eigenvalues] = eig(cov_matrix);
eigenvalues = diag(eigenvalues);

figure;
scatter3(EEG_Data(1, :), EEG_Data(2, :), EEG_Data(3, :), 20, t, 'filled');
colormap('winter');
colorbar;
hold on;

% Plot eigenvectors
origin = mean(EEG_Data, 2); % Compute the mean as the origin for the vectors
scale_factor = max(eigenvalues) * 2; % Adjust the scale factor for longer vectors

for i = 1:size(eigenvectors, 2)
    quiver3(origin(1), origin(2), origin(3), eigenvectors(1, i)*eigenvalues(1)*4, eigenvectors(2, i)*eigenvalues(2)*4, eigenvectors(3, i)*eigenvalues(3)*4, 'LineWidth', 2, 'Color', 'red');
end

hold off;

title('EEG Data with Eigenvectors');
xlabel('EEG Channel 1');
ylabel('EEG Channel 2');
zlabel('EEG Channel 3');
grid on;

% whitenning
whitening_transform = diag(1./sqrt(eigenvalues)) * eigenvectors';

% Apply whitening transform to EEG data
whitened_EEG_Data = whitening_transform * EEG_Data;
figure;
for channel = 1:size(whitened_EEG_Data, 1)
    subplot(size(EEG_Data, 1), 1, channel);
    plot(t, whitened_EEG_Data(channel, :));
    title(['Channel ' num2str(channel)]);
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    xlim([t(1) t(end)]); 
end
figure;
scatter3(whitened_EEG_Data(1, :), whitened_EEG_Data(2, :), whitened_EEG_Data(3, :), 20, t, 'filled');
colormap('winter');
colorbar;

title('Whitened EEG Data');
xlabel('Whitened EEG Channel 1');
ylabel('Whitened EEG Channel 2');
zlabel('Whitened EEG Channel 3');
grid on;

%% Part 4
[coeff, score, latent] = pca(EEG_Data');

whitening_transform = diag(1./sqrt(latent)) * coeff';

whitened_EEG_Data = whitening_transform * EEG_Data;

% Visualize whitened data in time domain
figure;
for channel = 1:size(whitened_EEG_Data, 1)
    subplot(size(EEG_Data, 1), 1, channel);
    plot(t, whitened_EEG_Data(channel, :));
    title(['Channel ' num2str(channel)]);
    xlabel('Time (seconds)');
    ylabel('Amplitude');
    xlim([t(1) t(end)]); 
end

% Visualize whitened data in 3D scatter plot
figure;
scatter3(whitened_EEG_Data(1, :), whitened_EEG_Data(2, :), whitened_EEG_Data(3, :), 20, t, 'filled');
colormap('winter');
colorbar;

title('Whitened EEG Data');
xlabel('Whitened EEG Channel 1');
ylabel('Whitened EEG Channel 2');
zlabel('Whitened EEG Channel 3');
grid on;

%% Part 5

% Perform Singular Value Decomposition (SVD)
[U, S, V] = svd(EEG_Data', 'econ');

% Calculate whitening transformation matrix
whitening_transform = U * diag(1./sqrt(diag(S))) * V';

% Apply whitening transform to EEG data
whitened_EEG_Data = whitening_transform * EEG_Data;

% Visualize elongations using RQ and LQ
figure;
subplot(1, 2, 1);
quiver3(zeros(1, size(V, 2)), zeros(1, size(V, 2)), zeros(1, size(V, 2)), V(1, :), V(2, :), V(3, :), 'LineWidth', 2);
title('Elongations using RQ');
xlabel('EEG Channel 1');
ylabel('EEG Channel 2');
zlabel('EEG Channel 3');
grid on;

subplot(1, 2, 2);
for i = 1:size(U, 2)
    quiver3(0, 0, 0, U(1, i), U(2, i), U(3, i), 'LineWidth', 2);
    hold on;
end
title('Elongations using LQ');
xlabel('EEG Channel 1');
ylabel('EEG Channel 2');
zlabel('EEG Channel 3');
grid on;
hold off;

EEG_Data_T = EEG_Data_C';
[U, S, V] = svd(EEG_Data_T, 'econ');
X_white = U * diag(1 ./ sqrt(diag(S))) * V';

figure;
scatter3(X_white(:, 1), X_white(:, 2), X_white(:, 3), 20, t,'filled');
colormap('winter');
colorbar;
title('Whitened Data Scatter Plot (3D)');
xlabel('Whitened Feature 1');
ylabel('Whitened Feature 2');
zlabel('Whitened Feature 3');
grid on;
