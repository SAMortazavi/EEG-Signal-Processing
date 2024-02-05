close all
clear
clc

%% Loading
Data=load("Q1.mat");
X_Org=Data.X_org;
X1=Data.X1;
X2=Data.X2;
X3=Data.X3;
X4=Data.X4;
freq=100;
%% Part1
P=400;
len=length(X_Org(1,:));
X_shifted= zeros([8,10000]);
X_shifted(:,P+1:len)=X_Org(:,1:len-P);
X_shifted(:,1:P)=X_Org(:,len-P+1:10000);
Px=cov(X_Org*X_shifted');
Px_hat=0.5*(Px + Px');
Cx=cov(X_Org');
[V_priodic,Lan_priodic]=eig(Px_hat,Cx);
[d,ind] = sort(diag(Lan_priodic),'descend');
w_S1=V_priodic(:,ind(1));
S_priodic=w_S1'*X_Org;
X_hat_S1= w_S1*S_priodic;
Error_X_S1=norm(X_hat_S1-X1);
offset = max(abs(X_hat_S1(:)));
disp_eeg(X_hat_S1,offset,freq);
offset = max(abs(X1(:)));
disp_eeg(X1,offset,freq);
%% Part 2 
P_test=300;
i=0;
Error_min=10000;
while i<400
    X_shifted_B= zeros([8,10000]);
    X_shifted_B(:,P_test+i+1:len)=X_Org(:,1:len-(P_test+i));
    X_shifted_B(:,1:P_test+i)=X_Org(:,len-(P_test+i)+1:10000);
    Px_B=cov(X_Org*X_shifted');
    Px_B_hat=0.5*(Px_B + Px_B');
    [V_priodic,Lan_priodic]=eig(Px_B_hat,Cx);
    [d,ind] = sort(diag(Lan_priodic),'descend');
    w_S1=V_priodic(:,ind(1));
    S_priodic_B=w_S1'*X_Org;
    X_hat_S1_B= w_S1*S_priodic_B;
    Error_X_S1_B=norm(X_hat_S1_B-X1);
    if Error_X_S1_B < Error_min
        P_min=P_test+i;
        Error_min=Error_X_S1_B;
        W_min=w_S1;
        S_priodic_min=S_priodic_B;
        X_hat_S1_min=X_hat_S1_B;
    end
    i=i+1;
end
offset = max(abs(X_hat_S1_min(:)));
disp_eeg(X_hat_S1_min,offset,freq);
offset = max(abs(X1(:)));
disp_eeg(X1,offset,freq);
%% Part 3
T_on= Data.T1;
X_on=zeros([8,10000]);
for row=1:8
    X_on(row,:)=X_Org(row,:).*T_on;
end
Cx_on=cov(X_on');
[V_nonstationary,Lan_nonstationar]=eig(Cx_on,Cx);
[d,ind] = sort(diag(Lan_nonstationar),'descend');
W_S2=V_nonstationary(:,ind(1));
S2_hat=W_S2'*X_Org;
X_S2_hat=W_S2*S2_hat;
Error_X_S2=norm(X_S2_hat-X2);

offset = max(abs(X_S2_hat(:)));
disp_eeg(X_S2_hat,offset,freq);
offset = max(abs(X2(:)));
disp_eeg(X2,offset,freq);
%% Part 4 
T_on_2=Data.T2;
X_on_2=zeros([8,10000]);
for row=1:8
    X_on_2(row,:)=X_Org(row,:).*T_on_2;
end
Cx_on_2=cov(X_on_2');
[V_nonstationary_2,Lan_nonstationar_2]=eig(Cx_on_2,Cx);
[d,ind] = sort(diag(Lan_nonstationar_2),'descend');
W_S2_2=V_nonstationary_2(:,ind(1));
S2_hat_2=W_S2_2'*X_Org;
X_S2_hat_2=W_S2_2*S2_hat_2;
Error_X_S2_2=norm(X_S2_hat_2-X2);

offset = max(abs(X_S2_hat_2(:)));
disp_eeg(X_S2_hat_2,offset,freq);
offset = max(abs(X2(:)));
disp_eeg(X2,offset,freq);
%% Part 5
Fs = 100;
Fpass = [10, 15];

filterOrder = 100; 

bpFilt = designfilt('bandpassfir', 'FilterOrder', filterOrder, 'CutoffFrequency1', Fpass(1), 'CutoffFrequency2', Fpass(2), 'SampleRate', Fs, 'Window', @kaiser);
% Apply the filter to the data
X_Filtered = filtfilt(bpFilt, X_Org');

X_Filtered = X_Filtered';

Px_Filtered=cov(X_Filtered*X_Filtered');
[V_Filtered,Lan_Filtered]=eig(Px_Filtered,Cx);
[d,ind] = sort(diag(Lan_Filtered),'descend');
W_S3=V_Filtered(:,ind(1));
S3_hat=W_S3'*X_Org;
X_S3_hat=W_S3*S3_hat;
Error_X_S3=norm(X_S3_hat-X3);

offset = max(abs(X_S3_hat(:)));
disp_eeg(X_S3_hat,offset,freq);
offset = max(abs(X3(:)));
disp_eeg(X3,offset,freq);

%% Part 6

Band_width=4;

Error_min=10000;

while Band_width < 15

    for i= 5:(25-Band_width)
        Fpass = [i, i+Band_width];
        bpFilt = designfilt('bandpassfir', 'FilterOrder', filterOrder, 'CutoffFrequency1', Fpass(1), 'CutoffFrequency2', Fpass(2), 'SampleRate', Fs, 'Window', @kaiser);
        X_Filtered_H = filtfilt(bpFilt, X_Org');
        X_Filtered_H = X_Filtered_H';
        Px_Filtered_H=cov(X_Filtered*X_Filtered_H');
        [V_Filtered_H,Lan_Filtered_H]=eig(Px_Filtered_H,Cx);
        [d,ind] = sort(diag(Lan_Filtered_H),'descend');

    end
    W_S3_H=V_Filtered_H(:,ind(1));
    S3_hat_H=W_S3_H'*X_Org;
    X_S3_hat_H=W_S3_H*S3_hat_H;
    Error_X_S3_H=norm(X_S3_hat_H-X3);
    if Error_X_S3_H < Error_min
        Error_min=Error_X_S3_H;
        Band_width_min=Band_width;
        range_Min=Fpass;
    end

    Band_width=Band_width+1;

end

offset = max(abs( X_S3_hat_H(:)));
disp_eeg( X_S3_hat_H,offset,freq);
offset = max(abs(X3(:)));
disp_eeg(X3,offset,freq);

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