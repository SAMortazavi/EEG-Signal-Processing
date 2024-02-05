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

%% Part 2  and 3 GEVD
%%  Noise 3 -10db

EEG_N310_On  =  T_on_rebuilt .* EEG_N3_10db;
Cx_hat_GEVD_N310  =  cov(EEG_N310_On');
Cx_N310  =  cov( EEG_N3_10db');
[V_GEVD_N310,Lan_GEVD_N310]  =  eig(Cx_hat_GEVD_N310 , Cx_N310);
[d,ind] = sort(diag(Lan_GEVD_N310),'descend');
W_GEVD_1_N310 = V_GEVD_N310(:,ind(1));
W_GEVD_2_N310=V_GEVD_N310(:,ind(2));
S1_GEVD_N310 = W_GEVD_1_N310'*EEG_N3_10db;
X_den_N310 =  W_GEVD_1_N310 * S1_GEVD_N310;
Error_X_GEVD_N310 = norm(X_den_N310-X_Org);



%Figures
% Noisy and Pure
Electrodes={'Pure 13','Pure 24','Noisy 13','Noisy 24'};
Showing_List=[X_Org(13,:);X_Org(24,:);EEG_N3_10db(13,:);EEG_N3_10db(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure and Noisy data N3 -10db');
% Denoised
Electrodes={'Denoised 13','Deniosed 24'};
Showing_List=[X_den_N310(13,:);X_den_N310(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Denoised data N3 -10dbdata N3 -10db');


%% Noise 3 -20db

EEG_N320_On  =  T_on_rebuilt .* EEG_N3_20db;
Cx_hat_GEVD_N320  =  cov( EEG_N320_On');
Cx_N320  =  cov( EEG_N3_20db');
[V_GEVD_N320,Lan_GEVD_N320]  =  eig(Cx_hat_GEVD_N320 , Cx_N320);
[d,ind] = sort(diag(Lan_GEVD_N320),'descend');
W_GEVD_1_N320 = V_GEVD_N320(:,ind(32));
W_GEVD_2_N320=V_GEVD_N320(:,ind(2));
S1_GEVD_N320 = W_GEVD_1_N320'*EEG_N3_20db;
X_den_N320 =  W_GEVD_1_N320 * S1_GEVD_N320;
Error_X_GEVD_N320 = norm(X_den_N320-X_Org);



%Figures
% Noisy and Pure
Electrodes={'Pure 13','Pure 24','Noisy 13','Noisy 24'};
Showing_List=[X_Org(13,:);X_Org(24,:);EEG_N3_20db(13,:);EEG_N3_20db(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure and Noisy data N3 -20db');
% Denoised
Electrodes={'Denoised 13','Deniosed 24'};
Showing_List=[X_den_N320(13,:);X_den_N320(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Denoised data N3 -10dbdata N3 -20db');


%% Noise 4 -10db

EEG_N410_On  =  T_on_rebuilt .* EEG_N4_10db;
Cx_hat_GEVD_N410  =  cov( EEG_N410_On');
Cx_N410  =  cov( EEG_N4_10db');
[V_GEVD_N410,Lan_GEVD_N410]  =  eig(Cx_hat_GEVD_N410 , Cx_N410);
[d,ind] = sort(diag(Lan_GEVD_N410),'descend');
W_GEVD_1_N410 = V_GEVD_N410(:,ind(1));
W_GEVD_2_N410=V_GEVD_N410(:,ind(2));
S1_GEVD_N410 = W_GEVD_1_N410'*EEG_N4_10db;
X_den_N410 =  W_GEVD_1_N410 * S1_GEVD_N410;
Error_X_GEVD_N410 = norm(X_den_N410-X_Org);



%Figures
% Noisy and Pure
Electrodes={'Pure 13','Pure 24','Noisy 13','Noisy 24'};
Showing_List=[X_Org(13,:);X_Org(24,:);EEG_N4_10db(13,:);EEG_N4_10db(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure and Noisy data N4 -10db');
% Denoised
Electrodes={'Denoised 13','Deniosed 24'};
Showing_List=[X_den_N410(13,:);X_den_N410(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Denoised data N3 -10dbdata N4 -10db');

%% Noise 4 -20db
EEG_N420_On  =  T_on_rebuilt .* EEG_N4_20db;
Cx_hat_GEVD_N420  =  cov( EEG_N420_On');
Cx_N420  =  cov( EEG_N4_20db');
[V_GEVD_N420,Lan_GEVD_N420]  =  eig(Cx_hat_GEVD_N420 , Cx_N420);
[d,ind] = sort(diag(Lan_GEVD_N420),'descend');
W_GEVD_1_N420 = V_GEVD_N420(:,ind(1));
S1_GEVD_N420 = W_GEVD_1_N420'*EEG_N4_20db;
X_den_N420 =  W_GEVD_1_N420 * S1_GEVD_N420;
Error_X_GEVD_N420 = norm(X_den_N420-X_Org);



%Figures
% Noisy and Pure
Electrodes={'Pure 13','Pure 24','Noisy 13','Noisy 24'};
Showing_List=[X_Org(13,:);X_Org(24,:);EEG_N4_20db(13,:);EEG_N4_20db(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure and Noisy data N4 -20db');
% Denoised
Electrodes={'Denoised 13','Deniosed 24'};
Showing_List=[X_den_N420(13,:);X_den_N420(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Denoised data N3 -20db');

%% Part 4
Numerator_3_10=0;
Numerator_3_20=0;
Numerator_4_10=0;
Numerator_4_20=0;
Denominator=0;
for i=1:32
    for j=1:10240
        Numerator_3_10=Numerator_3_10+(X_Org(i,j)-X_den_N310(i,j))^2;
        Numerator_3_20=Numerator_3_20+(X_Org(i,j)-X_den_N320(i,j))^2;
        Numerator_4_10=Numerator_4_10+(X_Org(i,j)-X_den_N410(i,j))^2;
        Numerator_4_20=Numerator_4_20+(X_Org(i,j)-X_den_N420(i,j))^2;
        Denominator=Denominator+(X_Org(i,j))^2;
    end
end
RRMSE_N3_10db=sqrt(Numerator_3_10/Denominator);
RRMSE_N3_20db=sqrt(Numerator_3_20/Denominator);
RRMSE_N4_10db=sqrt(Numerator_4_10/Denominator);
RRMSE_N4_20db=sqrt(Numerator_4_20/Denominator);






%% Functions
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