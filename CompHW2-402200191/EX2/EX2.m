close all
clear
clc;
%% part 1 
Data=load("Ex2.mat");
Data_EEG=Data.X_org;
Noise_3=Data.X_noise_3;
Noise_4=Data.X_noise_4;
Electrodes=load('Electrodes.mat').Electrodes.labels;
Power_Data_EEG = mean(abs(Data_EEG(:)).^2);
Power_Noise_3 = mean(abs(Noise_3(:)).^2);
Power_Noise_4 = mean(abs(Noise_4(:)).^2);
feq = 200 ;
% Noise 3
SNR10=-10;
sigma_Noise3_10=sqrt(Power_Data_EEG*10^(-SNR10/10)/(Power_Noise_3));
EEG_N3_10db= Data_EEG+sigma_Noise3_10*Noise_3;

SNR20=-20;
sigma_Noise3_20=sqrt(Power_Data_EEG*10^(-SNR20/10)/(Power_Noise_3));
EEG_N3_20db= Data_EEG+sigma_Noise3_20*Noise_3;


% Noise 3
SNR10=-10;
sigma_Noise4_10=sqrt(Power_Data_EEG*10^(-SNR10/10)/(Power_Noise_4));
EEG_N4_10db= Data_EEG+sigma_Noise4_10*Noise_4;

SNR20=-20;
sigma_Noise4_20=sqrt(Power_Data_EEG*10^(-SNR20/10)/(Power_Noise_4));
EEG_N4_20db= Data_EEG+sigma_Noise4_20*Noise_4;

% figures 
%Noise 3 10db
offset = max(abs(EEG_N3_10db(:))) ;
disp_eeg(EEG_N3_10db,offset,feq,Electrodes,'Noise 3 with -10db SNR');

%Noise 3 20db

offset = max(abs(EEG_N3_20db(:))) ;
disp_eeg(EEG_N3_20db,offset,feq,Electrodes,'Noise 3 with -20db SNR');

%Noise 4 10db

offset = max(abs(EEG_N4_10db(:))) ;
disp_eeg(EEG_N4_10db,offset,feq,Electrodes,'Noise 4 with -10db SNR');

%Noise 4 20db
offset = max(abs(EEG_N4_20db(:))) ;
disp_eeg(EEG_N4_20db,offset,feq,Electrodes,'Noise 4 with -20db SNR');



%% Part 2
[F,W,K]=COM2R(EEG_N3_10db,32);
source_coeff_N3_10=W;
obs_coeff_N3_10=F;
[F,W,K]=COM2R(EEG_N3_20db,32);
source_coeff_N3_20=W;
obs_coeff_N3_20=F;

[F,W,K]=COM2R(EEG_N4_10db,32);
source_coeff_N4_10=W;
obs_coeff_N4_10=F;
[F,W,K]=COM2R(EEG_N4_20db,32);
source_coeff_N4_20=W;
obs_coeff_N4_20=F;

% Sources

sources_N3_10db=source_coeff_N3_10*EEG_N3_10db;
sources_N3_20db=source_coeff_N3_20*EEG_N3_20db;
sources_N4_10db=source_coeff_N4_10*EEG_N4_10db;
sources_N4_20db=source_coeff_N4_20*EEG_N4_20db;


offset = max(abs(sources_N3_10db(:))) ;
disp_eeg(sources_N3_10db,offset,feq,Electrodes,'Sources For N3 -10db');

offset = max(abs(sources_N3_20db(:))) ;
disp_eeg(sources_N3_20db,offset,feq,Electrodes,'Sources For N3 -20db');

offset = max(abs(sources_N4_10db(:))) ;
disp_eeg(sources_N4_10db,offset,feq,Electrodes,'Sources For N4 -10db');

offset = max(abs(sources_N4_20db(:))) ;
disp_eeg(sources_N4_20db,offset,feq,Electrodes,'Sources For N4 -20db');

% PCA
[X_wh_N310,invWh_N310, Wh_N310]=MyPCA(EEG_N3_10db');
[X_wh_N320,invWh_N320, Wh_N320]=MyPCA(EEG_N3_20db');
[X_wh_N410,invWh_N410, Wh_N410]=MyPCA(EEG_N4_10db');
[X_wh_N420,invWh_N420, Wh_N420]=MyPCA(EEG_N4_20db');

offset = max(abs(X_wh_N310(:))) ;
disp_eeg(X_wh_N310,offset,feq,Electrodes,'Sources For N3 -10db with PCA');

offset = max(abs(X_wh_N320(:))) ;
disp_eeg(X_wh_N320,offset,feq,Electrodes,'Sources For N3 -20db with PCA');
offset = max(abs(X_wh_N410(:))) ;
disp_eeg(X_wh_N410,offset,feq,Electrodes,'Sources For N4 -10db with PCA');

offset = max(abs(X_wh_N420(:))) ;
disp_eeg(X_wh_N420,offset,feq,Electrodes,'Sources For N4 -20db with PCA');
%% Part 3,4
X_Den_N3_10= obs_coeff_N3_10(:,[4,12]) * sources_N3_10db([4,12],:);
X_Den_N3_20= obs_coeff_N3_20(:,[8,17,18]) * sources_N3_10db([8,17,18],:);
X_Den_N4_10= obs_coeff_N4_10(:,[2,3,27]) * sources_N4_10db([2,3,27],:);
X_Den_N4_20= obs_coeff_N4_20(:,25) * sources_N4_20db(25,:);

X_Den_N3_10_PCA= invWh_N310(:,[8,10,18,32]) * X_wh_N310([8,10,18,32],:);
X_Den_N3_20_PCA= invWh_N320(:,[8,18,32]) * X_wh_N320([8,18,32],:);
X_Den_N4_10_PCA= invWh_N410(:,[13,17,19,30]) * X_wh_N410([13,17,19,30],:);
X_Den_N4_20_PCA= invWh_N420(:,[17,19,30]) * X_wh_N420([17,19,30],:);

%% Part 5 


X_Den_N3_10_D= X_Den_N3_10;
X_Den_N3_20_D= X_Den_N3_20;
X_Den_N4_10_D= X_Den_N4_10;
X_Den_N4_20_D= X_Den_N3_10;
% Figures
Electrodes={'Pure 13','Pure 24','Noisy 13','Noisy 24','Denoised 13','Deniosed 24'};
Showing_List=[Data_EEG(13,:);Data_EEG(24,:);EEG_N3_10db(13,:);EEG_N3_10db(24,:);X_Den_N3_10_D(13,:);X_Den_N3_10_D(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure,Noisy and denoised data N3 -10db');

Showing_List=[Data_EEG(13,:);Data_EEG(24,:);EEG_N3_20db(13,:);EEG_N3_20db(24,:);X_Den_N3_20_D(13,:);X_Den_N3_20_D(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure,Noisy and denoised data N3 -20db');


Showing_List=[Data_EEG(13,:);Data_EEG(24,:);EEG_N4_10db(13,:);EEG_N4_10db(24,:);X_Den_N4_10_D(13,:);X_Den_N4_10_D(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure,Noisy and denoised data N4 -10 db');

Showing_List=[Data_EEG(13,:);Data_EEG(24,:);EEG_N4_20db(13,:);EEG_N4_20db(24,:);X_Den_N4_20_D(13,:);X_Den_N4_20_D(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure,Noisy and denoised data N4 -20 db');

% PCA
Showing_List=[Data_EEG(13,:);Data_EEG(24,:);EEG_N3_10db(13,:);EEG_N3_10db(24,:);X_Den_N3_10_PCA(13,:);X_Den_N3_10_PCA(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure,Noisy and denoised data N3 -10db with PCA');

Showing_List=[Data_EEG(13,:);Data_EEG(24,:);EEG_N3_20db(13,:);EEG_N3_20db(24,:);X_Den_N3_20_PCA(13,:);X_Den_N3_20_PCA(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure,Noisy and denoised data N3 -20db with PCA');


Showing_List=[Data_EEG(13,:);Data_EEG(24,:);EEG_N4_10db(13,:);EEG_N4_10db(24,:);X_Den_N4_10_PCA(13,:);X_Den_N4_10_PCA(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure,Noisy and denoised data N4 -10 db with PCA');

Showing_List=[Data_EEG(13,:);Data_EEG(24,:);EEG_N4_20db(13,:);EEG_N4_20db(24,:);X_Den_N4_20_PCA(13,:);X_Den_N4_20_PCA(24,:)];
offset = max(abs(Showing_List(:))) ;
disp_eeg(Showing_List,offset,feq,Electrodes,'Pure,Noisy and denoised data N4 -20 db with PCA');

%% Part 6
Numerator_3_10=0;
Numerator_3_20=0;
Numerator_4_10=0;
Numerator_4_20=0;
Denominator=0;
for i=1:32
    for j=1:10240
        Numerator_3_10=Numerator_3_10+(Data_EEG(i,j)-X_Den_N3_10_D(i,j))^2;
        Numerator_3_20=Numerator_3_20+(Data_EEG(i,j)-X_Den_N3_20_D(i,j))^2;
        Numerator_4_10=Numerator_4_10+(Data_EEG(i,j)-X_Den_N4_10_D(i,j))^2;
        Numerator_4_20=Numerator_4_20+(Data_EEG(i,j)-X_Den_N4_20_D(i,j))^2;
        Denominator=Denominator+(Data_EEG(i,j))^2;
    end
end
RRMSE_N3_10db=sqrt(Numerator_3_10/Denominator);
RRMSE_N3_20db=sqrt(Numerator_3_20/Denominator);
RRMSE_N4_10db=sqrt(Numerator_4_10/Denominator);
RRMSE_N4_20db=sqrt(Numerator_4_20/Denominator);

Numerator_3_10=0;
Numerator_3_20=0;
Numerator_4_10=0;
Numerator_4_20=0;
Denominator=0;
for i=1:32
    for j=1:10240
        Numerator_3_10=Numerator_3_10+(Data_EEG(i,j)-X_Den_N3_10_PCA(i,j))^2;
        Numerator_3_20=Numerator_3_20+(Data_EEG(i,j)-X_Den_N3_20_PCA(i,j))^2;
        Numerator_4_10=Numerator_4_10+(Data_EEG(i,j)-X_Den_N4_10_PCA(i,j))^2;
        Numerator_4_20=Numerator_4_20+(Data_EEG(i,j)-X_Den_N4_20_PCA(i,j))^2;
        Denominator=Denominator+(Data_EEG(i,j))^2;
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

function [F,W,K]=COM2R(Y,Pest)
disp('COM2')
% Comon, version 6 march 92
% English comments added in 1994
% [F,delta]=aci(Y)
% Y is the observations matrix
% This routine outputs a matrix F such that Y=F*Z, Z=pinv(F)*Y,
% and components of Z are uncorrelated and approximately independent
% F is Nxr, where r is the dimension Z;
% Entries of delta are sorted in decreasing order;
% Columns of F are of unit norm;
% The entry of largest modulus in each column of F is positive real.
% Initial and final values of contrast can be fprinted for checking.
% REFERENCE: P.Comon, "Independent Component Analysis, a new concept?",
% Signal Processing, Elsevier, vol.36, no 3, April 1994, 287-314.
%
[N,TT]=size(Y);T=max(N,TT);N=min(N,TT);
if TT==N, Y=Y';[N,T]=size(Y);end; % Y est maintenant NxT avec N<T.
%%%% STEPS 1 & 2: whitening and projection (PCA)
[U,S,V]=svd(Y',0);tol=max(size(S))*norm(S)*eps;
s=diag(S);I=find(s<tol);

%--- modif de Laurent le 03/02/2009
r = min(Pest,N);
U=U(:,1:r);
S=S(1:r,1:r);
V=V(:,1:r);
%---

Z=U'*sqrt(T);L=V*S'/sqrt(T);F=L; %%%%%% on a Y=L*Z;
%%%%%% INITIAL CONTRAST
T=length(Z);contraste=0;
for i=1:r,
 gii=Z(i,:)*Z(i,:)'/T;Z2i=Z(i,:).^2;;giiii=Z2i*Z2i'/T;
 qiiii=giiii/gii/gii-3;contraste=contraste+qiiii*qiiii;
end;
%%%% STEPS 3 & 4 & 5: Unitary transform
S=Z;
if N==2,K=1;else,K=1+round(sqrt(N));end;  % K= max number of sweeps
Rot=eye(r);
for k=1:K,                           %%%%%% strating sweeps
Q=eye(r);
  for i=1:r-1,
  for j= i+1:r,
    S1ij=[S(i,:);S(j,:)];
    [Sij,qij]=tfuni4(S1ij);    %%%%%% processing a pair
    S(i,:)=Sij(1,:);S(j,:)=Sij(2,:);
    Qij=eye(r);Qij(i,i)=qij(1,1);Qij(i,j)=qij(1,2);
    Qij(j,i)=qij(2,1);Qij(j,j)=qij(2,2);
    Q=Qij*Q;
  end;
  end;
Rot=Rot*Q';
end;                                    %%%%%% end sweeps
F=F*Rot;
%%%%%% FINAL CONTRAST
S=Rot'*Z;
T=length(S);contraste=0;
for i=1:r,
 gii=S(i,:)*S(i,:)'/T;S2i=S(i,:).^2;;giiii=S2i*S2i'/T;
 qiiii=giiii/gii/gii-3;contraste=contraste+qiiii*qiiii;
end;
%%%% STEP 6: Norming columns
delta=diag(sqrt(sum(F.*conj(F))));
%%%% STEP 7: Sorting
[d,I]=sort(-diag(delta));E=eye(r);P=E(:,I)';delta=P*delta*P';F=F*P';
%%%% STEP 8: Norming
F=F*inv(delta);
%%%% STEP 9: Phase of columns
[y,I]=max(abs(F));
for i=1:r,Lambda(i)=conj(F(I(i),i));end;Lambda=Lambda./abs(Lambda);
F=F*diag(Lambda);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCUL DE LA MATRICE DE FILTRAGE
%---------------------------------
W = pinv(F);


end
function [X_whiten,invMat, whMat] = MyPCA(X,epsilon)
  if ~exist('epsilon','var')
    epsilon = 0.0001;
  end
  mu=mean(X(:));
  X=X-mu;
  A = X'*X;
  [V,D,notused] = svd(A);
  whMat = sqrt(size(X,1)-1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
  Xwh = X*whMat; 
  X_whiten=Xwh';
  invMat = pinv(whMat);
end