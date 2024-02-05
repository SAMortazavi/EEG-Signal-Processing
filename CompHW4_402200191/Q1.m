clc
clear
close all
%%  Loading 
EEG_ERP=load('ERP_EEG.mat');
EEG=EEG_ERP.ERP_EEG;
EEG=EEG.';
%% section a
t=0:1/240:1-1/240;
for N=100:100:2500
    EEG_mean=mean(EEG(1:N,:));
    plot(t,EEG_mean)
    xlabel('time')
    title('mean of N sample')
    grid on
    hold on
    legendInfo{N/100} = ['N= ' num2str(N)]; 
    
end
legend(legendInfo)

%% section b
EEG_max=zeros(1,2550);
EEG_max(1,1)=max(abs(EEG(1,:)));
for N=2:2550
    EEG_mean=mean(EEG(1:N,:));
    EEG_max(1,N)=max(abs(EEG_mean));
  
end
figure()
plot(EEG_max)
xlabel('N')
title('max of mean of N sample')
grid on
%% section c
EEG_rms=zeros(1,2550);
EEG_mean_i_1=EEG(1,:);
for N=2:2550
    EEG_mean=mean(EEG(1:N,:));
    EEG_rms_i=rms(EEG_mean-EEG_mean_i_1);
    EEG_rms(1,N)=EEG_rms_i;
    EEG_mean_i_1=EEG_mean;
end
figure()
plot(EEG_rms)
xlabel('N')
grid on
title('difference between rms')

x=500:2550;
figure()
plot(x,EEG_rms(1,500:end))
xlabel('N')
grid on
%% section e
EEG_mean=mean(EEG(1:1000,:));
EEG_mean1=mean(EEG(1:2550,:));
EEG_mean2=mean(EEG(1:333,:));
EEG_mean3= mean(EEG((randi(length(EEG),1,1000)),:));
EEG_mean4= mean(EEG((randi(length(EEG),1,333)),:));
figure()
plot(EEG_mean)
hold on
plot(EEG_mean1)
hold on
plot(EEG_mean2)
hold on
plot(EEG_mean3)
hold on
plot(EEG_mean4)
hold on
xlabel('N')
grid on
legend('N=1000','N=2550','N=N_0 /3','N=N_0 random','N=N_0/3 random')


