clc
clear
close all
%% Loading
EEG_SSVEP=load('SSVEP_EEG.mat');
EEG=EEG_SSVEP.SSVEP_Signal;
Events=EEG_SSVEP.Events;
Event_samples=EEG_SSVEP.Event_samples;
Fs=250;
%% section a1
EEG_filtered = zeros(6,117917);
for i=1:6
    EEG_filtered(i,:) = bandpass( EEG(i,:) , [ 1 40] , Fs);
end
%% section a2
EEG_ss=zeros(15,2);
for i=1:15
    EEG_ss(i,1)=Event_samples(1,i);
     EEG_ss(i,2)=250*5+Event_samples(1,i);
end
%% section a3
figure()
for i=1:15
subplot(3,5,i)
%pwelch
for j=1:6
pxx = pwelch(EEG_filtered(j,EEG_ss(i,1):EEG_ss(i,2)));
%pwelch(EEG_filtered(j,EEG_ss(i,1):EEG_ss(i,2)))
% plot(pow2db(pxx))
plot(pxx)
hold on
end
title(i+" -th trial")
ylabel("pow/freq")
xlabel("frequency")
grid on
legend(subplot(3,5,1),{'Pz','Qz','P7','P8','O2','O1'})

end
%% section 2b
f=zeros(1,5);
y=zeros(1,5);
for i=1:5
   f(1,i)= Events(1,3*(i-1)+1);
end
t=5*Fs;
y1=generator(f(1,1),t);
y2=generator(f(1,2),t);
y3=generator(f(1,3),t);
y4=generator(f(1,4),t);
y5=generator(f(1,5),t);
rmax=zeros(1,15);
for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(:,EEG_ss(i,1):EEG_ss(i,2)-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(:,EEG_ss(i,1):EEG_ss(i,2)-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(:,EEG_ss(i,1):EEG_ss(i,2)-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(:,EEG_ss(i,1):EEG_ss(i,2)-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(:,EEG_ss(i,1):EEG_ss(i,2)-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc=length(m)/15
%% section 2c
f=zeros(1,5);
y=zeros(1,5);
for i=1:5
   f(1,i)= Events(1,3*(i-1)+1);
end
t=5*Fs;
y1=generator(f(1,1),t);
y2=generator(f(1,2),t);
y3=generator(f(1,3),t);
y4=generator(f(1,4),t);
y5=generator(f(1,5),t);
rmax=zeros(1,15);
for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,2)-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,2)-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,2)-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,2)-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,2)-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_Pz=length(m)/15

for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(3,EEG_ss(i,1):EEG_ss(i,2)-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(3,EEG_ss(i,1):EEG_ss(i,2)-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(3,EEG_ss(i,1):EEG_ss(i,2)-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(3,EEG_ss(i,1):EEG_ss(i,2)-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(3,EEG_ss(i,1):EEG_ss(i,2)-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_P7=length(m)/15

for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(5,EEG_ss(i,1):EEG_ss(i,2)-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(5,EEG_ss(i,1):EEG_ss(i,2)-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(5,EEG_ss(i,1):EEG_ss(i,2)-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(5,EEG_ss(i,1):EEG_ss(i,2)-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(5,EEG_ss(i,1):EEG_ss(i,2)-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_O2=length(m)/15
%% section 2d

f=zeros(1,5);
y=zeros(1,5);
for i=1:5
   f(1,i)= Events(1,3*(i-1)+1);
end
t=3*Fs;
y1=generator(f(1,1),t);
y2=generator(f(1,2),t);
y3=generator(f(1,3),t);
y4=generator(f(1,4),t);
y5=generator(f(1,5),t);
rmax=zeros(1,15);
for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+3*Fs-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+3*Fs-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+3*Fs-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+3*Fs-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+3*Fs-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_3T=length(m)/15

t=2*Fs;
y1=generator(f(1,1),t);
y2=generator(f(1,2),t);
y3=generator(f(1,3),t);
y4=generator(f(1,4),t);
y5=generator(f(1,5),t);
rmax=zeros(1,15);
for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+2*Fs-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+2*Fs-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+2*Fs-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+2*Fs-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+2*Fs-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_2T=length(m)/15

t=Fs;
y1=generator(f(1,1),t);
y2=generator(f(1,2),t);
y3=generator(f(1,3),t);
y4=generator(f(1,4),t);
y5=generator(f(1,5),t);
rmax=zeros(1,15);
for i=1:15
    [~,~,r1] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+Fs-1).',y1.');
    [~,~,r2] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+Fs-1).',y2.');
    [~,~,r3] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+Fs-1).',y3.');
    [~,~,r4] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+Fs-1).',y4.');
    [~,~,r5] = canoncorr(EEG_filtered(1,EEG_ss(i,1):EEG_ss(i,1)+Fs-1).',y5.');
    r=[max(r1);max(r2);max(r3);max(r4);max(r5)];
    [n,~]=find(max(r)==r);
    rmax(i)=f(1,n);
end
[n,m]=find(Events==rmax);
acc_channel_T=length(m)/15

%% Function
function y=generator(f,t)
    num=floor(40/f);
    y=zeros(2*num,t);
    n=t/250;
    T=[0:1/250:n-1/250];
    for i=1:num
       y(2*i-1,:)= sin(2*pi*i*f*T);
       y(2*i,:)= cos(2*pi*i*f*T);
    end
end
