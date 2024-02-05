close all
clear 
clc
%% Part 1
Data = load("Ex2.mat");

X_Org = Data.X_org;

Noise_3 = Data.X_noise_3;
Noise_4 = Data.X_noise_4;

Power_X_Org  =  mean(abs(X_Org(:)).^2);
Power_Noise_3  =  mean(abs(Noise_3(:)).^2);
Power_Noise_4  =  mean(abs(Noise_4(:)).^2);
feq  =  200 ;
% Noise 3
SNR10 = -10;
sigma_Noise3_10 = sqrt(Power_X_Org*10^(-SNR10/10)/(Power_Noise_3));
EEG_N3_10db =  X_Org+sigma_Noise3_10*Noise_3;

SNR20 = -20;
sigma_Noise3_20 = sqrt(Power_X_Org*10^(-SNR20/10)/(Power_Noise_3));
EEG_N3_20db =  X_Org+sigma_Noise3_20*Noise_3;

SNR20 = -20;
sigma_Noise3_20 = sqrt(Power_X_Org*10^(-SNR20/10)/(Power_Noise_3));
EEG_N3_20db =  X_Org+sigma_Noise3_20*Noise_3;


% Noise 3
SNR10 = -10;
sigma_Noise4_10 = sqrt(Power_X_Org*10^(-SNR10/10)/(Power_Noise_4));
EEG_N4_10db =  X_Org+sigma_Noise4_10*Noise_4;

SNR20 = -20;
sigma_Noise4_20 = sqrt(Power_X_Org*10^(-SNR20/10)/(Power_Noise_4));
EEG_N4_20db =  X_Org+sigma_Noise4_20*Noise_4;
% Finding T_On
threshold = 0.59;
spike_indices = find(abs(X_Org) > threshold);
T_on = zeros(size(X_Org));
T_on(spike_indices) = 1;

column_sum = sum(T_on, 1);

row_threshold = 1/32* size(T_on, 1); 

spike_time_columns = find(column_sum > row_threshold);

T_on_rebuilt = zeros(1, size(X_Org, 2));
T_on_rebuilt(spike_time_columns) = 1;

%% Part 2 
%% Noise 3 -10db
W_N310=rand(32,1);
EEG_N3_10db_whitened= whiten(EEG_N3_10db');
for i=1:10000
    W_previous=W_N310;
    r=W_N310' * EEG_N3_10db_whitened;
    r_Plus=T_on_rebuilt .* r;
    for j=1:32
        W_N310(j)=sum(EEG_N3_10db_whitened(j,:).*r_Plus);
    end
    W_N310=W_N310/norm(W_N310);
    if sum(abs(W_N310)-abs(W_previous))==0
        break
    end
end
S_N310 = W_N310' * EEG_N3_10db_whitened;
X_Den_N310= W_N310 * S_N310;
Error_N310=norm(X_Den_N310 - X_Org);
%% Noise 3 -20
W_N320=rand(32,1);
EEG_N3_20db_whitened= whiten(EEG_N3_20db');
for i=1:500000
    W_previous=W_N320;
    r=W_N320' * EEG_N3_20db_whitened;
    r_Plus=T_on_rebuilt .* r;
    for j=1:32
        W_N320(j)=sum(EEG_N3_20db_whitened(j,:).*r_Plus);
    end
    W_N320=W_N320/norm(W_N320);
    if sum(abs(W_N320)-abs(W_previous))==0
        break
    end
end
S_N320 = W_N320' * EEG_N3_20db_whitened;
X_Den_N320 = W_N320 *S_N320;
Error_N320 = norm(X_Den_N320 - X_Org);

%% Noise 4 -10db
W_N410=rand(32,1);
EEG_N4_10db_whitened= whiten(EEG_N4_10db');
for i=1:10000
    W_previous=W_N410;
    r=W_N410' * EEG_N4_10db_whitened;
    r_Plus=T_on_rebuilt .* r;
    for j=1:32
        W_N410(j)=sum(EEG_N4_10db_whitened(j,:).*r_Plus);
    end
    W_N410=W_N410/norm(W_N410);
    if sum(abs(W_N410)-abs(W_previous))==0
        break
    end
end
S_N410 = W_N410' * EEG_N4_10db_whitened;
X_Den_N410= W_N410 * S_N410;
Error_N410=norm(X_Den_N410 - X_Org);

%% Noise 4 -20db
W_N420=rand(32,1);
EEG_N4_20db_whitened= whiten(EEG_N4_20db');
for i=1:500000
    W_previous=W_N420;
    r=W_N420' * EEG_N4_20db_whitened;
    r_Plus=T_on_rebuilt .* r;
    for j=1:32
        W_N420(j)=sum(EEG_N4_20db_whitened(j,:).*r_Plus);
    end
    W_N420=W_N420/norm(W_N420);
    if sum(abs(W_N420)-abs(W_previous))==0
        break
    end
end
S_N420 = W_N420' * EEG_N4_20db_whitened;
X_Den_N420= W_N420 * S_N420;
Error_N420=norm(X_Den_N420 - X_Org);

%% Part 3
%% N 310
Electrodes={'Pure 13','Pure 24','Noisy 13','Noisy 24'};
Showing_List=[X_Org(13,:);X_Org(24,:);EEG_N4_20db(13,:);EEG_N4_20db(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure and Noisy data N3 -10db');
% Denoised
Electrodes={'Denoised 13','Deniosed 24'};
Showing_List=[X_Den_N310(13,:);X_Den_N310(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Denoised data N3 -10db');

%% N 320
Electrodes={'Pure 13','Pure 24','Noisy 13','Noisy 24'};
Showing_List=[X_Org(13,:);X_Org(24,:);EEG_N4_20db(13,:);EEG_N4_20db(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure and Noisy data N3 -20db');
% Denoised
Electrodes={'Denoised 13','Deniosed 24'};
Showing_List=[X_Den_N320(13,:);X_Den_N320(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Denoised data N3 -20db');
%% N 410
Electrodes={'Pure 13','Pure 24','Noisy 13','Noisy 24'};
Showing_List=[X_Org(13,:);X_Org(24,:);EEG_N4_20db(13,:);EEG_N4_20db(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure and Noisy data N4 -10db');
% Denoised
Electrodes={'Denoised 13','Deniosed 24'};
Showing_List=[X_Den_N410(13,:);X_Den_N410(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Denoised data N4 -10db');

%% N 420
Electrodes={'Pure 13','Pure 24','Noisy 13','Noisy 24'};
Showing_List=[X_Org(13,:);X_Org(24,:);EEG_N4_20db(13,:);EEG_N4_20db(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure and Noisy data N4 -20db');
% Denoised
Electrodes={'Denoised 13','Deniosed 24'};
Showing_List=[X_Den_N420(13,:);X_Den_N420(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Denoised data N4 -20db');


%% Part 4
Numerator_3_10=0;
Numerator_3_20=0;
Numerator_4_10=0;
Numerator_4_20=0;
Denominator=0;
for i=1:32
    for j=1:10240
        Numerator_3_10=Numerator_3_10+(X_Org(i,j)-X_Den_N310(i,j))^2;
        Numerator_3_20=Numerator_3_20+(X_Org(i,j)-X_Den_N320(i,j))^2;
        Numerator_4_10=Numerator_4_10+(X_Org(i,j)-X_Den_N410(i,j))^2;
        Numerator_4_20=Numerator_4_20+(X_Org(i,j)-X_Den_N420(i,j))^2;
        Denominator=Denominator+(X_Org(i,j))^2;
    end
end
RRMSE_N3_10db=sqrt(Numerator_3_10/Denominator);
RRMSE_N3_20db=sqrt(Numerator_3_20/Denominator);
RRMSE_N4_10db=sqrt(Numerator_4_10/Denominator);
RRMSE_N4_20db=sqrt(Numerator_4_20/Denominator);



%% Functions
function [Xwh, mu, invMat, whMat] = whiten(X,epsilon)

if ~exist('epsilon','var')
    epsilon = 0.0001;
end
mu = mean(X); 
X = bsxfun(@minus, X, mu);
A = X'*X;
[V,D,notused] = svd(A);
whMat = sqrt(size(X,1)-1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
Xwh = X*whMat; 
Xwh=Xwh';
invMat = pinv(whMat);
end

function t = disp_eeg(X,offset,feq,ElecName,titre)
% function t = disp_eeg(X,offset,feq,ElecName,titre)
%
% inputs
%     X: dynamics to display. (nbchannels x nbsamples) matrix
%     offset: offset between channels (default max(abs(X)))
%     feq: sapling frequency (default 1)
%     ElecName: cell array of electrode labels (default {S1,S2,...})
%     titre: title of the figure
%
% output
%     t: time vector
%
% G. Birot 2010-02


%% Check arguments
[N K] = size(X);

if nargin < 4
    for n = 1:N
        ElecName{n}  = ['S',num2str(n)];
    end
    titre = [];
end

if nargin < 5
    titre = [];
end

if isempty(feq)
    feq = 1;
end

if isempty(ElecName)
    for n = 1:N
        ElecName{n}  = ['S',num2str(n)];
    end
end

if isempty(offset)
    offset = max(abs(X(:)));
end


%% Build dynamic matrix with offset and time vector
X = X + repmat(offset*(0:-1:-(N-1))',1,K);
t = (1:K)/feq;
graduations = offset*(0:-1:-(N-1))';
shiftvec = N:-1:1;
Ysup = max(X(1,:)) + offset;
Yinf = min(X(end,:)) - offset;
% YLabels = cell(N+2) ElecName(shiftvec)

%% Display
figure1 = figure;
% a1 = axes('YAxisLocation','right');
a2 = axes('YTickLabel',ElecName(shiftvec),'YTick',graduations(shiftvec),'FontSize',7);
ylim([Yinf Ysup]);
box('on');
grid('on')
hold('all');
plot(t,X');
xlabel('Time (seconds)','FontSize',10);
ylabel('Channels','FontSize',10);
title(titre);
hold off

end