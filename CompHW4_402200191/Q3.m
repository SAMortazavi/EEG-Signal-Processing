clc
clear
close all
%% Loading
EEG=load('Ex3.mat');
TestData=EEG.TestData;
TrainData=EEG.TrainData;
TrainLabel=EEG.TrainLabel;

%% CSP
n=find(TrainLabel);
Cx1=0;
for i=1:n
Cx1=Cx1+TrainData(:,:,n(i))*TrainData(:,:,n(i)).'/trace(TrainData(:,:,n(i))*TrainData(:,:,n(i)).');
end
Cx1=Cx1/length(n);

m=find(TrainLabel==0);
Cx2=0;
for i=1:m
Cx2=Cx2+TrainData(:,:,m(i))*TrainData(:,:,m(i)).'./trace(TrainData(:,:,m(i))*TrainData(:,:,m(i)).');
end
Cx2=Cx2/length(m);

% GEVD
[U,D]=eig(Cx1,Cx2);
D = diag(D);
[D, index] = sort(D, 'descend');
U = (U(:,index));
W=[U(:,1) U(:,end)];
class0=W.'*TrainData(:,:,m(1,1));

class1=W.'*TrainData(:,:,n(1,1));
figure()
subplot(2,1,1)
plot(class0(1,:))
hold on
plot(class1(1,:))
title('First Filter')
legend('class0','class1')
grid on

subplot(2,1,2)
plot(class0(2,:))
hold on
plot(class1(2,:))
title('Last Filter')
legend('class0','class1')
grid on
%% section b
Electrodes=load('AllElectrodes.mat');
elecX=zeros(30,1);
elecY=zeros(30,1);
label=[];

i=[37 7 5 38 40 42 10 47 45  15 13 48 50 52 18 32 55 23 22 21 20 31 57  58 59 60 26 25 27 64];
for j=1:30
    elecX(j)=Electrodes.AllElectrodes(i(j)).X;
    elecY(j)=Electrodes.AllElectrodes(i(j)).Y;
    label{j}=Electrodes.AllElectrodes(i(j)).labels;
    
end
subplot(1,2,1)
plottopomap(elecX,elecY, label,abs(U(:,1)))
title('First Filter')
subplot(1,2,2)
plottopomap(elecX,elecY, label,abs(U(:,end)))
title('Last Filter')

%% section 3 

a=[1 42 83 124 165];

for filter_num=1:15
    acc=0;
    for i=1:3
    test_4fold=TrainData(:,:,a(i):a(i+1));
    test_4fold_label=TrainLabel(1,a(i):a(i+1));
    train_4fold=TrainData;
    train_4fold(:,:,a(i):a(i+1))=[];
    train_4fold_label=TrainLabel;
    train_4fold_label(:,a(i):a(i+1))=[];
    n=find(train_4fold_label);
    Cx1=0;
    for j=1:n
        Cx1=Cx1+train_4fold(:,:,n(j))*train_4fold(:,:,n(j)).'/trace(train_4fold(:,:,n(j))*train_4fold(:,:,n(j)).');
    end
    Cx1=Cx1/length(n);

    m=find(train_4fold_label==0);
    Cx2=0;
    for j=1:m
        Cx2=Cx2+train_4fold(:,:,m(j))*train_4fold(:,:,m(j)).'./trace(train_4fold(:,:,m(j))*train_4fold(:,:,m(j)).');
    end
    Cx2=Cx2/length(m);

    % GEVD
    [U,D]=eig(Cx1,Cx2);
    D = diag(D);
    [D, index] = sort(D, 'descend');
    U = (U(:,index));
    W=[U(:,1:filter_num) U(:,end-filter_num+1:end)];
    
    var_class0=zeros(length(m),2*filter_num);
    for j=1: length(m)
        class0=W.'*train_4fold(:,:,m(1,j));
        var_class0(j,:)=var(class0');
    end
    
    var_class1=zeros(length(n),2*filter_num);
    for j=1: length(n)
        class1=W.'*train_4fold(:,:,n(1,j));
        var_class1(j,:)=var(class1');
    end
    
    var_test=zeros(length(test_4fold(1,1,:)),2*filter_num);
    for j=1: length(test_4fold(1,1,:))
        class_test=W.'*test_4fold(:,:,j);
        var_test(j,:)=var(class_test');
    end
    
    feature=[var_class0 ; var_class1];
    label=[zeros(length(m),1); ones(length(n),1)];
    Mdl = fitcknn(feature,label);
    label_test_predict = predict(Mdl,var_test);
    

    acc0=length(find(label_test_predict'==test_4fold_label))/length(test_4fold_label);
    acc=acc+acc0;
    end
  accuracy(filter_num)=acc/3;
end
x=1:15;
plot(2*x,accuracy)
xlabel('number of filter')
ylabel('acc')
grid on
%% section d
%best filter num
filter_num=6;
acc=0;
TestData=EEG.TestData;
TrainData=EEG.TrainData;
TrainLabel=EEG.TrainLabel;
    n=find(TrainLabel);
    Cx1=0;
    for j=1:n
        Cx1=Cx1+TrainData(:,:,n(j))*TrainData(:,:,n(j)).'/trace(TrainData(:,:,n(j))*TrainData(:,:,n(j)).');
    end
    Cx1=Cx1/length(n);

    m=find(TrainLabel==0);
    Cx2=0;
    for j=1:m
        Cx2=Cx2+TrainData(:,:,m(j))*TrainData(:,:,m(j)).'./trace(TrainData(:,:,m(j))*TrainData(:,:,m(j)).');
    end
    Cx2=Cx2/length(m);

    % GEVD
    [U,D]=eig(Cx1,Cx2);
    D = diag(D);
    [D, index] = sort(D, 'descend');
    U = (U(:,index));
    W=[U(:,1:filter_num) U(:,end-filter_num+1:end)];
    
    var_class0=zeros(length(m),2*filter_num);
    for j=1: length(m)
        class0=W.'*TrainData(:,:,m(1,j));
        var_class0(j,:)=var(class0');
    end
    
    var_class1=zeros(length(n),2*filter_num);
    for j=1: length(n)
        class1=W.'*TrainData(:,:,n(1,j));
        var_class1(j,:)=var(class1');
    end
    
    var_test=zeros(length(TestData(1,1,:)),2*filter_num);
    for j=1: length(TestData(1,1,:))
        class_test=W.'*TestData(:,:,j);
        var_test(j,:)=var(class_test');
    end
    
    feature=[var_class0 ; var_class1];
    label=[zeros(length(m),1); ones(length(n),1)];
    Mdl = fitcknn(feature,label);
    label_test_predict = predict(Mdl,var_test);
 %% save   
writematrix(label_test_predict,'Best_predict_Lables.xlsx');

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