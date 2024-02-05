close all
clear 
clc
%% Part 1
New_Data1=load("NewData1.mat").EEG_Sig;
New_Data2=load("NewData2.mat").EEG_Sig;
New_Data3=load("NewData3.mat").EEG_Sig;
New_Data4=load("NewData4.mat").EEG_Sig;
ElectrodesX=load("Electrodes.mat").Electrodes.X;
ElectrodesY=load("Electrodes.mat").Electrodes.Y;
ElectrodesZ=load("Electrodes.mat").Electrodes.Z;
Electrodes=load("Electrodes.mat").Electrodes.labels;
freq=250;

offset = max(abs(New_Data1(:))) ;
disp_eeg(New_Data1,offset,freq,Electrodes,'New Data 1');

offset = max(abs(New_Data2(:))) ;
disp_eeg(New_Data2,offset,freq,Electrodes,'New Data 2');

offset = max(abs(New_Data3(:))) ;
disp_eeg(New_Data3,offset,freq,Electrodes,'New Data 3');

offset = max(abs(New_Data4(:))) ;
disp_eeg(New_Data4,offset,freq,Electrodes,'New Data 4');



%% Part 3
[F,W,K]=COM2R(New_Data1,21);
source_coeff_D1=W;
obs_coeff_D1=F;

[F,W,K]=COM2R(New_Data2,21);
source_coeff_D2=W;
obs_coeff_D2=F;

[F,W,K]=COM2R(New_Data3,21);
source_coeff_D3=W;
obs_coeff_D3=F;


[F,W,K]=COM2R(New_Data4,21);
source_coeff_D4=W;
obs_coeff_D4=F;

Sources_ND1=source_coeff_D1*New_Data1;
Sources_ND2=source_coeff_D2*New_Data2;
Sources_ND3=source_coeff_D3*New_Data3;
Sources_ND4=source_coeff_D4*New_Data4;


%% Part 4
offset = max(abs(Sources_ND1(:))) ;
disp_eeg(Sources_ND1,offset,freq,Electrodes,'Sources ND1');

offset = max(abs(Sources_ND2(:))) ;
disp_eeg(Sources_ND2,offset,freq,Electrodes,'Sources ND2');

offset = max(abs(Sources_ND3(:))) ;
disp_eeg(Sources_ND3,offset,freq,Electrodes,'Sources ND3');

offset = max(abs(Sources_ND4(:))) ;
disp_eeg(Sources_ND4,offset,freq,Electrodes,'Sources ND4');


% Set parameters for pwelch
fs = 250;
window = hamming(256); 
noverlap = 128;

figure;

for i = 1:4
    subplot(2, 2, i);
    
    psd_values = zeros(size(eval(['Sources_ND', num2str(i)]), 1), numel(window)/2 + 1);

    for j = 1:size(eval(['Sources_ND', num2str(i)]), 1)
        [pxx, f] = pwelch(eval(['Sources_ND', num2str(i), '(j, :)']), window, noverlap, [], fs);
        psd_values(j, :) = pxx';
    end

    plot(f, mean(psd_values, 1));
    
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    title(['Mean PSD for Sources ND', num2str(i)]);
end

%% Special
% I wrote a function to show data for every columns.
plottopomap_all_sources(ElectrodesX, ElectrodesY, Electrodes, obs_coeff_D1);
plottopomap_all_sources(ElectrodesX, ElectrodesY, Electrodes, obs_coeff_D2);
plottopomap_all_sources(ElectrodesX, ElectrodesY, Electrodes, obs_coeff_D3);
plottopomap_all_sources(ElectrodesX, ElectrodesY, Electrodes, obs_coeff_D4);


%% Part 6
Selsource_ND1=[2,5,6,8,9,12];
Selsource_ND2=[1,6,16,19,21];
Selsource_ND3=[2,3,6,9,20,21];
Selsource_ND4=[6,9,20,21];

X_denoised_N1=obs_coeff_D1(:,Selsource_ND1)*Sources_ND1(Selsource_ND1,:);
X_denoised_N2=obs_coeff_D2(:,Selsource_ND2)*Sources_ND2(Selsource_ND2,:);
X_denoised_N3=obs_coeff_D3(:,Selsource_ND3)*Sources_ND3(Selsource_ND3,:);
X_denoised_N4=obs_coeff_D4(:,Selsource_ND4)*Sources_ND4(Selsource_ND4,:);

offset = max(abs(X_denoised_N1(:))) ;
disp_eeg(X_denoised_N1,offset,freq,Electrodes,'X denoised N1');

offset = max(abs(X_denoised_N2(:))) ;
disp_eeg(X_denoised_N2,offset,freq,Electrodes,'X denoised N2');

offset = max(abs(X_denoised_N3(:))) ;
disp_eeg(X_denoised_N3,offset,freq,Electrodes,'X denoised N3');

offset = max(abs(X_denoised_N4(:))) ;
disp_eeg(X_denoised_N4,offset,freq,Electrodes,'X denoised N4');



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


function plottopomap(elocsX,elocsY,elabels,data)

% define XY points for interpolation
interp_detail = 100;
interpX = linspace(min(elocsX)-.2,max(elocsX)+.25,interp_detail);
interpY = linspace(min(elocsY),max(elocsY),interp_detail);

% meshgrid is a function that creates 2D grid locations based on 1D inputs
[gridX,gridY] = meshgrid(interpX,interpY);
% Interpolate the data on a 2D grid
interpFunction = TriScatteredInterp(elocsY,elocsX,data);
topodata = interpFunction(gridX,gridY);

% plot map
contourf(interpY,interpX,topodata);
hold on
scatter(elocsY,elocsX,10,'ro','filled');
for i=1:length(elocsX)
    text(elocsY(i),elocsX(i),elabels(i))
end
set(gca,'xtick',[])
set(gca,'ytick',[])
end


function plottopomap_all_sources(ElectrodesX, ElectrodesY, Electrodes, obs_coeff_ND)
    num_of_col=length(obs_coeff_ND(1,:));
    figure;
    for i=1:num_of_col
    subplot(7,3,i)
    plottopomap(ElectrodesX,ElectrodesY,Electrodes,obs_coeff_ND(:,i))

    end

end
