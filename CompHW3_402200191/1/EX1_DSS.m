clear
close all
clc
%% Loading
Data=load("Q1.mat");
X_Org=Data.X_org;
X1=Data.X1;
X2=Data.X2;
X3=Data.X3;
X4=Data.X4;
freq=100;
%% Part 1
X_whitened=whiten(X_Org');
P=400;
L=length(X_whitened(1,:))/P;
Num_S=8;
w_S1=rand(Num_S,1);
i=0;
while i<10000
    w_previous=w_S1;
    r_S1=w_S1'*X_whitened;
    r_S1_hat=zeros([1,P]);
    for j=1:L
        r_S1_hat=r_S1_hat+r_S1(1+(j-1)*P:P*j);
    end
    r_S1_hat=r_S1_hat/P;
    r_S1_plus=zeros([1,length(X_whitened(1,:))]);
    for j=1:L
        r_S1_plus(1+(j-1)*P:P*j)=r_S1_hat;
    end
    for j=1:Num_S
        w_S1(j)=sum(X_whitened(j,:).*r_S1_plus);
    end
    w_S1=w_S1/norm(w_S1);
    if sum(abs(w_S1)-abs(w_previous))==0
        break
    end
    i=i+1;
end
S_priodic=w_S1'*X_whitened;
X_hat_S1=w_S1*S_priodic;
Error_X_S1=norm(X_hat_S1-X1);

offset = max(abs(X_hat_S1(:)));
disp_eeg(X_hat_S1,offset,freq);

offset = max(abs(X1(:)));
disp_eeg(X1,offset,freq);
%% Part 2
Periods=[400,500,625];
k=1;
Error_min=10000;
while k<=3
    P=Periods(k);
    P
    L_B=length(X_whitened(1,:))/P;
    w_B=rand(Num_S,1);
    i=0;
    while i<10000
        w_previous_B=w_B;
        r_S1=w_B'*X_whitened;
        r_S1_hat=zeros([1,P]);
        for j=1:L_B
            r_S1_hat=r_S1_hat+r_S1(1+(j-1)*P:P*j);
        end
        r_S1_hat=r_S1_hat/P;
        r_S1_plus=zeros([1,length(X_whitened(1,:))]);
        for j=1:L_B
            r_S1_plus(1+(j-1)*P:P*j)=r_S1_hat;
        end
        for j=1:Num_S
            w_B(j)=sum(X_whitened(j,:).*r_S1_plus);
        end
        w_B=w_B/norm(w_B);
        if sum(abs(w_B)-abs(w_previous_B))==0
            break
        end
        i=i+1;
    end
    S_priodic_B=w_B'*X_whitened;
    X_hat_S1_B=w_B*S_priodic_B;
    Error_X_S1_B=norm(X_hat_S1_B-X1)
    if Error_X_S1_B<Error_min
        S_priodic_min=S_priodic_B;
        P_min=P;
        X_hat_S1_min=X_hat_S1_B;
        Error_min=Error_X_S1_B;
    end
    k=k+1;
end

offset = max(abs(X_hat_S1_min(:)));
disp_eeg(X_hat_S1_min,offset,freq);

offset = max(abs(X1(:)));
disp_eeg(X1,offset,freq);
%% Part 3
W_S2=rand(Num_S,1);
T_on=Data.T1;
i=0;
while i<10000
    w_previous_C=W_S2;
    r_C=W_S2'*X_whitened;
    r_S2_plus=r_C.*T_on;
    for j=1:Num_S
        W_S2(j)=sum(X_whitened(j,:).*r_S2_plus);
    end
    W_S2=W_S2/norm(W_S2);
    if sum(abs(W_S2)-abs(w_previous_C))==0
        break
    end


   i=i+1;
end

S_nonstationary=W_S2'*X_whitened;
X_hat_S2=W_S2*S_nonstationary;
Error_X_S2=norm(X_hat_S2-X2);


offset = max(abs(X_hat_S2(:)));
disp_eeg(X_hat_S2,offset,freq);

offset = max(abs(X2(:)));
disp_eeg(X2,offset,freq);

%% Part 4
W_S2_2=rand(Num_S,1);
T_on_2=Data.T2;
for i=1:10000
    w_previous_D=W_S2_2;
    r_D=W_S2_2'*X_whitened;
    r_S2_plus_D=r_D.*T_on_2;  
    for j=1:Num_S
        W_S2_2(j)=sum(X_whitened(j,:).*r_S2_plus_D);
    end
    W_S2_2=W_S2_2/norm(W_S2_2);
    if sum(abs(W_S2_2)-abs(w_previous_D))==0
        break
    end
end
S_nonstationary_2=W_S2_2'*X_whitened;
X_hat_S2_2=W_S2_2*S_nonstationary_2;
Error_X_S2_2=norm(X_hat_S2_2-X2);

offset = max(abs(X_hat_S2_2(:)));
disp_eeg(X_hat_S2_2,offset,freq);

offset = max(abs(X2(:)));
disp_eeg(X2,offset,freq);

%% Part 5
Fs = 100;
Fpass = [10, 15];
% Design a bandpass filter
filterOrder = 100; % Adjust the filter order based on your requirements

bpFilt = designfilt('bandpassfir', 'FilterOrder', filterOrder, 'CutoffFrequency1', Fpass(1), 'CutoffFrequency2', Fpass(2), 'SampleRate', Fs, 'Window', @kaiser);

W_S3=rand(Num_S,1);
for i=1 : 10000
    w_previous_H=W_S3;
    r_S3= W_S3' *X_whitened;
    r_S3_filtered = filtfilt(bpFilt, r_S3');
    r_S3_Plus = r_S3_filtered';
    for j=1:Num_S
        W_S3(j)=sum(X_whitened(j,:).*r_S3_Plus);
    end
      W_S3=W_S3/norm(W_S3);
     if sum(abs(W_S3)-abs(w_previous_H))==0
        break
    end
end
S_3=W_S3'*X_whitened;
X_hat_S3=W_S3*S_3;
Error_X_S3=norm(X_hat_S3-X3);


offset = max(abs(X_hat_S3(:)));
disp_eeg(X_hat_S3,offset,freq);

offset = max(abs(X3(:)));
disp_eeg(X3,offset,freq);

%% Part 6

Band_width=4;

Error_min_H=10000;

while Band_width < 15
        
    for k= 5:(25-Band_width)
        Fpass = [k, k+Band_width];
        bpFilt = designfilt('bandpassfir', 'FilterOrder', filterOrder, 'CutoffFrequency1', Fpass(1), 'CutoffFrequency2', Fpass(2), 'SampleRate', Fs, 'Window', @kaiser);
        W_S3_H=rand(Num_S,1);
        for i=1 : 10000
            w_previous_H=W_S3_H;
            r_S3_H= W_S3_H' *X_whitened;
            r_S3_filtered_H = filtfilt(bpFilt, r_S3_H');
            r_S3_Plus_H = r_S3_filtered_H';
            for j=1:Num_S
                W_S3_H(j)=sum(X_whitened(j,:).*r_S3_Plus_H);
            end
            W_S3_H=W_S3_H/norm(W_S3_H);
            if sum(abs(W_S3_H)-abs(w_previous_H))==0
                break
            end
        end
        S_3_H=W_S3_H'*X_whitened;
        X_hat_S3_H=W_S3_H*S_3_H;
        Error_X_S3_H=norm(X_hat_S3_H-X3);
        if Error_X_S3_H < Error_min_H
            Error_min_H=Error_X_S3_H;
            Band_width_min=Band_width;
            range_Min=Fpass;
        end
        k
    end



    Band_width=Band_width+1;

end
offset = max(abs(X_hat_S3_H(:)));
disp_eeg(X_hat_S3_H,offset,freq);

offset = max(abs(X3(:)));
disp_eeg(X3,offset,freq);

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

%% Function
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