clear all
close all
clc
%% section a
load('ElecPosXYZ.mat') ;
%Forward Matrix
ModelParams.R = [8 8.5 9.2] ;
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3];
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];

Resolution = 1 ;
[LocMat,GainMat] = ForwardModel_3shell(Resolution, ModelParams) ;
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'r','.');
xlabel('x')
ylabel('y')
zlabel('z')
title('Location of dipoles')
%% section b
Elec_name = {};
Elec_Pos = zeros(3, 21);

for i = 1:21
    Elec_Pos(:, i) = ElecPos{1, i}.XYZ;
    Elec_name{i} = ElecPos{1, i}.Name;
end

scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'r', '.');
xlabel('x')
ylabel('y')
zlabel('z')
title('Location of dipoles with electrodes')
hold on
scatter3(ModelParams.R(3) * Elec_Pos(1, :), ModelParams.R(3) * Elec_Pos(2, :), ModelParams.R(3) * Elec_Pos(3, :), 'b', '.');
for i = 1:21
    text(ModelParams.R(3) * Elec_Pos(1, i), ModelParams.R(3) * Elec_Pos(2, i), ModelParams.R(3) * Elec_Pos(3, i), Elec_name{i});
end

%% section c
clc;
x = -5; y = 0; z = 6.2;
dipole_index = 1200;

vector_normalize = LocMat(:, dipole_index) / sqrt(sum(LocMat(:, dipole_index).^2, 'all'));

subplot(1, 2, 1)
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'r', '.');
hold on
scatter3(LocMat(1, dipole_index), LocMat(2, dipole_index), LocMat(3, dipole_index), 'g', '*');
hold on
plot3([vector_normalize(1, 1) + LocMat(1, dipole_index), LocMat(1, dipole_index)], ...
      [vector_normalize(2, 1) + LocMat(2, dipole_index), LocMat(2, dipole_index)], ...
      [vector_normalize(3, 1) + LocMat(3, dipole_index), LocMat(3, dipole_index)], '->')
xlabel('x')
ylabel('y')
zlabel('z')
title('Location of dipoles with electrodes only with e_q')
hold on
scatter3(ModelParams.R(3) * Elec_Pos(1, :), ModelParams.R(3) * Elec_Pos(2, :), ModelParams.R(3) * Elec_Pos(3, :), 'b', '.');

for i = 1:21
    text(ModelParams.R(3) * Elec_Pos(1, i), ModelParams.R(3) * Elec_Pos(2, i), ModelParams.R(3) * Elec_Pos(3, i), Elec_name{i});
end

subplot(1, 2, 2)
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'r', '.');
hold on
scatter3(LocMat(1, dipole_index), LocMat(2, dipole_index), LocMat(3, dipole_index), 'g', '*');
hold on
plot3([ModelParams.R(1) * vector_normalize(1, 1) + LocMat(1, dipole_index), LocMat(1, dipole_index)], ...
      [ModelParams.R(1) * vector_normalize(2, 1) + LocMat(2, dipole_index), LocMat(2, dipole_index)], ...
      [ModelParams.R(1) * vector_normalize(3, 1) + LocMat(3, dipole_index), LocMat(3, dipole_index)], '->')
xlabel('x')
ylabel('y')
zlabel('z')
title('Location of dipoles with electrodes with q*e_q')
hold on
scatter3(ModelParams.R(3) * Elec_Pos(1, :), ModelParams.R(3) * Elec_Pos(2, :), ModelParams.R(3) * Elec_Pos(3, :), 'b', '.');

for i = 1:21
    text(ModelParams.R(3) * Elec_Pos(1, i), ModelParams.R(3) * Elec_Pos(2, i), ModelParams.R(3) * Elec_Pos(3, i), Elec_name{i});
end

%% section d
Interictal=load('Matlab/Interictal.mat') ;
Interictal=Interictal.Interictal;
M=(GainMat(:,3*(dipole_index-1)+1:3*dipole_index)*vector_normalize)*Interictal(10,:);

figure()
for i=1:21
    subplot( 3 ,7 ,i)
    plot( M(i,:))
    xlabel('time')
    title( Elec_name{i})
end
%% section e
clc;
index_peak=zeros(21,8);
peak=zeros(21,8);
for i=1:21
    index_peak(i,1)=find(abs(M(i,1:1500))==max(abs(M(i,1:1500))));
    peak(i,1)=M(i,index_peak(i,1));
    
    index_peak(i,2)=1500+find(abs(M(i,1500:5000))==max(abs(M(i,1500:5000))));
    peak(i,2)=M(i,index_peak(i,2));
    
    index_peak(i,3)=5000+find(abs(M(i,5000:5500))==max(abs(M(i,5000:5500))));
    peak(i,3)=M(i,index_peak(i,3));
    
    index_peak(i,4)=5500+find(abs(M(i,5500:7000))==max(abs(M(i,5500:7000))));
    peak(i,4)=M(i,index_peak(i,4));
    
    index_peak(i,5)=7000+find(abs(M(i,7000:8000))==max(abs(M(i,7000:8000))));
    peak(i,5)=M(i,index_peak(i,5));
    
    index_peak(i,6)=8000+find(abs(M(i,8000:9000))==max(abs(M(i,8000:9000))));
    peak(i,6)=M(i,index_peak(i,6));
    
    index_peak(i,7)=9000+find(abs(M(i,9000:10000))==max(abs(M(i,9000:10000))));
    peak(i,7)=M(i,index_peak(i,7));
    
    index_peak(i,8)=10000+find(abs(M(i,10000:end))==max(abs(M(i,10000:end))));
    peak(i,8)=M(i,index_peak(i,8));
end
figure()
for i=1:21
    subplot( 3 ,7 ,i)
    plot( M(i,:))
    xlabel('time')
    title( Elec_name{i})
    hold on
    scatter( index_peak(i,:), peak(i,:),'r','.')
    
end

mean_peak=zeros(21,1);
index=[index_peak(:,1)-3:index_peak(:,1)+3 index_peak(:,2)-3:index_peak(:,2)+3 index_peak(:,3)-3:index_peak(:,3)+3 index_peak(:,4)-3:index_peak(:,4)+3 index_peak(:,5)-3:index_peak(:,5)+3 index_peak(:,6)-3:index_peak(:,6)+3 index_peak(:,7)-3:index_peak(:,7)+3 index_peak(:,8)-3:index_peak(:,8)+3];
for i=1:21
    mean_peak(i,1)=mean(abs(M(i,index_peak(1,:))),'all');
end
figure()
scatter3(ModelParams.R(3)*Elec_Pos(1,:),ModelParams.R(3)*Elec_Pos(2,:),ModelParams.R(3)*Elec_Pos(3,:),'*','CData',mean_peak);
hold on
colorbar
for i = 1:21
    text(ModelParams.R(3) * Elec_Pos(1, i), ModelParams.R(3) * Elec_Pos(2, i), ModelParams.R(3) * Elec_Pos(3, i), Elec_name{i});
end
%% section f
figure()
Display_Potential_3D(ModelParams.R(3),mean_peak)
%% section g
% MNE
alpha=0.1;
Q_MNE=GainMat.'*inv(GainMat*GainMat.'+alpha*eye(21))*mean_peak;
Q_MNE_hat=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    Q_MNE_hat(1,i)=sqrt(sum(Q_MNE((i-1)*3+1:i*3).^2,'all'));
end
index_dipole_MNE=find(Q_MNE_hat==max(Q_MNE_hat))
location_dipole_MNE=LocMat(:,index_dipole_MNE)
i=index_dipole_MNE;
vector_normalize_MNE=Q_MNE((i-1)*3+1:i*3)./sqrt(sum(Q_MNE((i-1)*3+1:i*3).^2,'all'))
%% WMNE
G=GainMat;
w=zeros(length(Q_MNE)/3,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
   w(i,i)=sqrt(sum( G(:,(i-1)*3+1:i*3)*G(:,(i-1)*3+1:i*3).','all'));
end
W=kron(w,eye(3));
Q_WMNE=inv(W.'*W)*G.'*inv(G*(inv(W.'*W))*G.'+alpha*eye(21))*mean_peak;
Q_WMNE_hat=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    Q_WMNE_hat(1,i)=sqrt(sum(Q_WMNE((i-1)*3+1:i*3).^2,'all'));
end
index_dipole_WMNE=find(Q_WMNE_hat==max(Q_WMNE_hat))
location_dipole_WMNE=LocMat(:,index_dipole_WMNE)
i=index_dipole_WMNE;
vector_normalize_WMNE=Q_WMNE((i-1)*3+1:i*3)./sqrt(sum(Q_WMNE((i-1)*3+1:i*3).^2,'all'))

%% loreta 
d=1;
A1=zeros(length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    for j=1:length(Q_MNE)/3
        a=sqrt(sum((LocMat(:,i)-LocMat(:,j)).^2,'all'));
        if a==d
            A1(i,j)=1/6;
        end
    end
end
A0=inv(diag(A1*ones(length(Q_MNE)/3,1)))*A1;
A=kron(A0,eye(3));
B=(6/d^2)*(A-eye(length(Q_MNE)));

G=GainMat;
w=zeros(length(Q_MNE)/3,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
   w(i,i)=sqrt(sum( G(:,(i-1)*3+1:i*3)*G(:,(i-1)*3+1:i*3).','all'));
end
W=kron(w,eye(3));
W_loreta=W*B.'*B*W;

Q_Loreta=inv(W_loreta.'*W_loreta)*G.'*inv(G*(inv(W_loreta.'*W_loreta))*G.'+alpha*eye(21))*mean_peak;
Q_Loreta_hat=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    Q_Loreta_hat(1,i)=sqrt(sum(Q_Loreta((i-1)*3+1:i*3).^2,'all'));
end
index_dipole_Loreta=find(Q_Loreta_hat==max(Q_Loreta_hat))
location_dipole_Loreta=LocMat(:,index_dipole_Loreta)
i=index_dipole_Loreta;
vector_normalize_Loreta=Q_Loreta((i-1)*3+1:i*3)./sqrt(sum(Q_Loreta((i-1)*3+1:i*3).^2,'all'))

%% section h
clc;
error_location_MNE=sqrt(sum((LocMat(:,dipole_index)-location_dipole_MNE).^2,'all'))
dot_vector_MNE=dot(vector_normalize,vector_normalize_MNE);
theta_vector_MNE_radian=acos(dot_vector_MNE/sqrt(sum((vector_normalize).^2,'all'))*sqrt(sum((vector_normalize_MNE).^2,'all')))
theta_vector_MNE_degree=rad2deg( theta_vector_MNE_radian ) 

error_location_WMNE=sqrt(sum((LocMat(:,dipole_index)-location_dipole_WMNE).^2,'all'))
dot_vector_WMNE=dot(vector_normalize,vector_normalize_WMNE);
theta_vector_WMNE_radian=acos(dot_vector_WMNE/sqrt(sum((vector_normalize).^2,'all'))*sqrt(sum((vector_normalize_WMNE).^2,'all')))
theta_vector_WMNE_degree=rad2deg( theta_vector_WMNE_radian ) 

error_location_Loreta=sqrt(sum((LocMat(:,dipole_index)-location_dipole_Loreta).^2,'all'))
dot_vector_Loreta=dot(vector_normalize,vector_normalize_Loreta);
theta_vector_Loreta_radian=acos(dot_vector_Loreta/sqrt(sum((vector_normalize).^2,'all'))*sqrt(sum((vector_normalize_Loreta).^2,'all')))
theta_vector_Loreta_degree=rad2deg( theta_vector_Loreta_radian ) 

%% section i
% deep dipole

%% section i-c
clc;
x=2 ;y=-1;z=1.2;
dipole_index_deep=512;
vector_normalize=LocMat(:,dipole_index_deep)./sqrt(sum(LocMat(:,dipole_index_deep).^2,'all'));
subplot(1,2,1)
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'r','.');
hold on
scatter3(LocMat(1,dipole_index_deep),LocMat(2,dipole_index_deep),LocMat(3,dipole_index_deep),'g','*');
hold on
plot3( [ vector_normalize(1,1)+LocMat(1,dipole_index_deep),LocMat(1,dipole_index_deep) ] ,[ vector_normalize(2,1)+LocMat(2,dipole_index_deep),LocMat(2,dipole_index_deep)] , [vector_normalize(3,1)+LocMat(3,dipole_index_deep),LocMat(3,dipole_index_deep)] ,'->')
xlabel('x')
ylabel('y')
zlabel('z')
title('Location of dipoles with electrodes only with e_q')
hold on
scatter3(ModelParams.R(3)*Elec_Pos(1,:),ModelParams.R(3)*Elec_Pos(2,:),ModelParams.R(3)*Elec_Pos(3,:),'b','.');
hold on
for i = 1:21
    text(ModelParams.R(3) * Elec_Pos(1, i), ModelParams.R(3) * Elec_Pos(2, i), ModelParams.R(3) * Elec_Pos(3, i), Elec_name{i});
end
subplot(1,2,2)
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'r','.');
hold on
scatter3(LocMat(1,dipole_index_deep),LocMat(2,dipole_index_deep),LocMat(3,dipole_index_deep),'g','*');
hold on
plot3( [ LocMat(1,dipole_index_deep) ,ModelParams.R(1)*vector_normalize(1,1)+LocMat(1,dipole_index_deep)] ,[LocMat(2,dipole_index_deep), ModelParams.R(1)*vector_normalize(2,1)+LocMat(2,dipole_index_deep)] , [LocMat(3,dipole_index_deep),ModelParams.R(1)*vector_normalize(3,1)+LocMat(3,dipole_index_deep)],'->')
xlabel('x')
ylabel('y')
zlabel('z')
title('Location of dipoles with electrodes with q*e_q')
hold on
scatter3(ModelParams.R(3)*Elec_Pos(1,:),ModelParams.R(3)*Elec_Pos(2,:),ModelParams.R(3)*Elec_Pos(3,:),'b','.');
hold on
for i = 1:21
    text(ModelParams.R(3) * Elec_Pos(1, i), ModelParams.R(3) * Elec_Pos(2, i), ModelParams.R(3) * Elec_Pos(3, i), Elec_name{i});
end
%% section i-d
Interictal=load('Matlab/Interictal.mat') ;
Interictal=Interictal.Interictal;
M=(GainMat(:,3*(dipole_index_deep-1)+1:3*dipole_index_deep)*vector_normalize)*Interictal(10,:);

figure()
for i=1:21
    subplot( 3 ,7 ,i)
    plot( M(i,:))
    xlabel('time')
    title( Elec_name{i})
end
%% section i-e
clc;
index_peak=zeros(21,8);
peak=zeros(21,8);
for i=1:21
    index_peak(i,1)=find(abs(M(i,1:1500))==max(abs(M(i,1:1500))));
    peak(i,1)=M(i,index_peak(i,1));
    
    index_peak(i,2)=1500+find(abs(M(i,1500:5000))==max(abs(M(i,1500:5000))));
    peak(i,2)=M(i,index_peak(i,2));
    
    index_peak(i,3)=5000+find(abs(M(i,5000:5500))==max(abs(M(i,5000:5500))));
    peak(i,3)=M(i,index_peak(i,3));
    
    index_peak(i,4)=5500+find(abs(M(i,5500:7000))==max(abs(M(i,5500:7000))));
    peak(i,4)=M(i,index_peak(i,4));
    
    index_peak(i,5)=7000+find(abs(M(i,7000:8000))==max(abs(M(i,7000:8000))));
    peak(i,5)=M(i,index_peak(i,5));
    
    index_peak(i,6)=8000+find(abs(M(i,8000:9000))==max(abs(M(i,8000:9000))));
    peak(i,6)=M(i,index_peak(i,6));
    
    index_peak(i,7)=9000+find(abs(M(i,9000:10000))==max(abs(M(i,9000:10000))));
    peak(i,7)=M(i,index_peak(i,7));
    
    index_peak(i,8)=10000+find(abs(M(i,10000:end))==max(abs(M(i,10000:end))));
    peak(i,8)=M(i,index_peak(i,8));
end
figure()
for i=1:21
    subplot( 3 ,7 ,i)
    plot( M(i,:))
    xlabel('time')
    title( Elec_name{i})
    hold on
    scatter( index_peak(i,:), peak(i,:),'r','.')
    
end

mean_peak=zeros(21,1);
index=[index_peak(:,1)-3:index_peak(:,1)+3 index_peak(:,2)-3:index_peak(:,2)+3 index_peak(:,3)-3:index_peak(:,3)+3 index_peak(:,4)-3:index_peak(:,4)+3 index_peak(:,5)-3:index_peak(:,5)+3 index_peak(:,6)-3:index_peak(:,6)+3 index_peak(:,7)-3:index_peak(:,7)+3 index_peak(:,8)-3:index_peak(:,8)+3];
for i=1:21
    mean_peak(i,1)=mean(abs(M(i,index_peak(1,:))),'all');
end
figure()
scatter3(ModelParams.R(3)*Elec_Pos(1,:),ModelParams.R(3)*Elec_Pos(2,:),ModelParams.R(3)*Elec_Pos(3,:),'*','CData',mean_peak);
hold on
colorbar
for i = 1:21
    text(ModelParams.R(3) * Elec_Pos(1, i), ModelParams.R(3) * Elec_Pos(2, i), ModelParams.R(3) * Elec_Pos(3, i), Elec_name{i});
end
%% section i-f
figure()
Display_Potential_3D(ModelParams.R(3),mean_peak)
%% section i-g
% MNE
alpha=0.1;
Q_MNE=GainMat.'*inv(GainMat*GainMat.'+alpha*eye(21))*mean_peak;
Q_MNE_hat=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    Q_MNE_hat(1,i)=sqrt(sum(Q_MNE((i-1)*3+1:i*3).^2,'all'));
end
index_dipole_MNE=find(Q_MNE_hat==max(Q_MNE_hat))
location_dipole_MNE=LocMat(:,index_dipole_MNE)
i=index_dipole_MNE;
vector_normalize_MNE=Q_MNE((i-1)*3+1:i*3)./sqrt(sum(Q_MNE((i-1)*3+1:i*3).^2,'all'))
%% WMNE
G=GainMat;
w=zeros(length(Q_MNE)/3,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
   w(i,i)=sqrt(sum( G(:,(i-1)*3+1:i*3)*G(:,(i-1)*3+1:i*3).','all'));
end
W=kron(w,eye(3));
Q_WMNE=inv(W.'*W)*G.'*inv(G*(inv(W.'*W))*G.'+alpha*eye(21))*mean_peak;
Q_WMNE_hat=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    Q_WMNE_hat(1,i)=sqrt(sum(Q_WMNE((i-1)*3+1:i*3).^2,'all'));
end
index_dipole_WMNE=find(Q_WMNE_hat==max(Q_WMNE_hat))
location_dipole_WMNE=LocMat(:,index_dipole_WMNE)
i=index_dipole_WMNE;
vector_normalize_WMNE=Q_WMNE((i-1)*3+1:i*3)./sqrt(sum(Q_WMNE((i-1)*3+1:i*3).^2,'all'))

%% section j
% loreta 
d=1;
A1=zeros(length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    for j=1:length(Q_MNE)/3
        a=sqrt(sum((LocMat(:,i)-LocMat(:,j)).^2,'all'));
        if a==d
            A1(i,j)=1/6;
        end
    end
end
A0=inv(diag(A1*ones(length(Q_MNE)/3,1)))*A1;
A=kron(A0,eye(3));
B=(6/d^2)*(A-eye(length(Q_MNE)));

G=GainMat;
w=zeros(length(Q_MNE)/3,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
   w(i,i)=sqrt(sum( G(:,(i-1)*3+1:i*3)*G(:,(i-1)*3+1:i*3).','all'));
end
W=kron(w,eye(3));
W_loreta=W*B.'*B*W;

Q_Loreta=inv(W_loreta.'*W_loreta)*G.'*inv(G*(inv(W_loreta.'*W_loreta))*G.'+alpha*eye(21))*mean_peak;
Q_Loreta_hat=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    Q_Loreta_hat(1,i)=sqrt(sum(Q_Loreta((i-1)*3+1:i*3).^2,'all'));
end
index_dipole_Loreta=find(Q_Loreta_hat==max(Q_Loreta_hat))
location_dipole_Loreta=LocMat(:,index_dipole_Loreta)
i=index_dipole_Loreta;
vector_normalize_Loreta=Q_Loreta((i-1)*3+1:i*3)./sqrt(sum(Q_Loreta((i-1)*3+1:i*3).^2,'all'))
%% section i-h
clc;
error_location_MNE=sqrt(sum((LocMat(:,dipole_index_deep)-location_dipole_MNE).^2,'all'))
dot_vector_MNE=dot(vector_normalize,vector_normalize_MNE);
theta_vector_MNE_radian=acos(dot_vector_MNE/sqrt(sum((vector_normalize).^2,'all'))*sqrt(sum((vector_normalize_MNE).^2,'all')))
theta_vector_MNE_degree=rad2deg( theta_vector_MNE_radian ) 

error_location_WMNE=sqrt(sum((LocMat(:,dipole_index_deep)-location_dipole_WMNE).^2,'all'))
dot_vector_WMNE=dot(vector_normalize,vector_normalize_WMNE);
theta_vector_WMNE_radian=acos(dot_vector_WMNE/sqrt(sum((vector_normalize).^2,'all'))*sqrt(sum((vector_normalize_WMNE).^2,'all')))
theta_vector_WMNE_degree=rad2deg( theta_vector_WMNE_radian ) 

error_location_Loreta=sqrt(sum((LocMat(:,dipole_index_deep)-location_dipole_Loreta).^2,'all'))
dot_vector_Loreta=dot(vector_normalize,vector_normalize_Loreta);
theta_vector_Loreta_radian=acos(dot_vector_Loreta/sqrt(sum((vector_normalize).^2,'all'))*sqrt(sum((vector_normalize_Loreta).^2,'all')))
theta_vector_Loreta_degree=rad2deg( theta_vector_Loreta_radian ) 

%% section m
clc;
index_dipoles=[65:69 80:84 95:99 110:114];
vector_normalize=LocMat(:,index_dipoles)./sqrt(sum(LocMat(:,index_dipoles).^2,'all'));
subplot(1,2,1)
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'r','.');
hold on
scatter3(LocMat(1,index_dipoles),LocMat(2,index_dipoles),LocMat(3,index_dipoles),'g','*');
hold on
for i=1:20
    plot3( [LocMat(1,index_dipoles(i)) , vector_normalize(1,i)+LocMat(1,index_dipoles(i))] ,[LocMat(2,index_dipoles(i)) , vector_normalize(2,i)+LocMat(2,index_dipoles(i))] , [LocMat(3,index_dipoles(i)) , vector_normalize(3,i)+LocMat(3,index_dipoles(i))] ,'->')
    %plot3( [ vector_normalize(1,i)+LocMat(1,dipole_index_deep),LocMat(1,dipole_index_deep) ] ,[ vector_normalize(2,1)+LocMat(2,dipole_index_deep),LocMat(2,dipole_index_deep)] , [vector_normalize(3,1)+LocMat(3,dipole_index_deep),LocMat(3,dipole_index_deep)] ,'->')

    hold on
end
xlabel('x')
ylabel('y')
zlabel('z')
title('Location of dipoles with electrodes only with e_q')
hold on
scatter3(ModelParams.R(3)*Elec_Pos(1,:),ModelParams.R(3)*Elec_Pos(2,:),ModelParams.R(3)*Elec_Pos(3,:),'b','.');
hold on
for i = 1:21
    text(ModelParams.R(3) * Elec_Pos(1, i), ModelParams.R(3) * Elec_Pos(2, i), ModelParams.R(3) * Elec_Pos(3, i), Elec_name{i});
end
subplot(1,2,2)
scatter3(LocMat(1,:),LocMat(2,:),LocMat(3,:),'r','.');
hold on
scatter3(LocMat(1,index_dipoles),LocMat(2,index_dipoles),LocMat(3,index_dipoles),'g','*');
hold on
for i=1:20
    plot3( [LocMat(1,index_dipoles(i)) , ModelParams.R(1)*vector_normalize(1,i)+LocMat(1,index_dipoles(i))] ,[LocMat(2,index_dipoles(i)) , ModelParams.R(1)*vector_normalize(2,i)+LocMat(2,index_dipoles(i))] , [LocMat(3,index_dipoles(i)) , ModelParams.R(1)*vector_normalize(3,i)+LocMat(3,index_dipoles(i))] ,'->')
    %plot3( [ vector_normalize(1,i)+LocMat(1,dipole_index_deep),LocMat(1,dipole_index_deep) ] ,[ vector_normalize(2,1)+LocMat(2,dipole_index_deep),LocMat(2,dipole_index_deep)] , [vector_normalize(3,1)+LocMat(3,dipole_index_deep),LocMat(3,dipole_index_deep)] ,'->')

    hold on
end
xlabel('x')
ylabel('y')
zlabel('z')
title('Location of dipoles with electrodes with q*e_q')
hold on
scatter3(ModelParams.R(3)*Elec_Pos(1,:),ModelParams.R(3)*Elec_Pos(2,:),ModelParams.R(3)*Elec_Pos(3,:),'b','.');
hold on
for i = 1:21
    text(ModelParams.R(3) * Elec_Pos(1, i), ModelParams.R(3) * Elec_Pos(2, i), ModelParams.R(3) * Elec_Pos(3, i), Elec_name{i});
end



%% section n- d
Interictal=load('Matlab/Interictal.mat') ;
Interictal=Interictal.Interictal;
M=(GainMat(:,3*(index_dipoles-1)+1:3*index_dipoles)*vector_normalize)*Interictal;
figure()
for i=1:21
    subplot( 3 ,7 ,i)
    plot( M(i,:))
    xlabel('time')
    title( Elec_name{i})
end
%% section n-e
clc;
index_peak=zeros(21,4);
peak=zeros(21,4);
for i=1:21
    index_peak(i,1)=find(abs(M(i,1:2000))==max(abs(M(i,1:2000))));
    peak(i,1)=M(i,index_peak(i,1));
    
    index_peak(i,2)=2000+find(abs(M(i,2000:6000))==max(abs(M(i,2000:6000))));
    peak(i,2)=M(i,index_peak(i,2));
    
    index_peak(i,3)=6000+find(abs(M(i,6000:9000))==max(abs(M(i,6000:9000))));
    peak(i,3)=M(i,index_peak(i,3));
    
    index_peak(i,4)=9000+find(abs(M(i,9000:end))==max(abs(M(i,9000:end))));
    peak(i,4)=M(i,index_peak(i,4));
end
figure()
for i=1:21
    subplot( 3 ,7 ,i)
    plot( M(i,:))
    xlabel('time')
    title( Elec_name{i})
    hold on
    scatter( index_peak(i,:), peak(i,:),'r','.')
    
end

mean_peak=zeros(21,1);
index=[index_peak(:,1)-3:index_peak(:,1)+3 index_peak(:,2)-3:index_peak(:,2)+3 index_peak(:,3)-3:index_peak(:,3)+3 index_peak(:,4)-3:index_peak(:,4)+3 ];
for i=1:21
    mean_peak(i,1)=mean(abs(M(i,index_peak(1,:))),'all');
end
figure()
scatter3(ModelParams.R(3)*Elec_Pos(1,:),ModelParams.R(3)*Elec_Pos(2,:),ModelParams.R(3)*Elec_Pos(3,:),'*','CData',mean_peak);
hold on
colorbar
for i = 1:21
    text(ModelParams.R(3) * Elec_Pos(1, i), ModelParams.R(3) * Elec_Pos(2, i), ModelParams.R(3) * Elec_Pos(3, i), Elec_name{i});
end
%% section n-f
figure()
Display_Potential_3D(ModelParams.R(3),mean_peak)

%% section L
% MNE
alpha=0.1;
Q_MNE=GainMat.'*inv(GainMat*GainMat.'+alpha*eye(21))*mean_peak;
Q_MNE_hat=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    Q_MNE_hat(1,i)=sqrt(sum(Q_MNE((i-1)*3+1:i*3).^2,'all'));
end
[Q index]=sort( Q_MNE_hat ,'ascend');
%% WMNE
G=GainMat;
w=zeros(length(Q_MNE)/3,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
   w(i,i)=sqrt(sum( G(:,(i-1)*3+1:i*3)*G(:,(i-1)*3+1:i*3).','all'));
end
W=kron(w,eye(3));
Q_WMNE=inv(W.'*W)*G.'*inv(G*(inv(W.'*W))*G.'+alpha*eye(21))*mean_peak;
Q_WMNE_hat=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    Q_WMNE_hat(1,i)=sqrt(sum(Q_WMNE((i-1)*3+1:i*3).^2,'all'));
end
[Q index]=sort( Q_WMNE_hat ,'ascend');
%% section o
% MNE
labels = zeros(1,length(Q_MNE)/3);
labels(index_dipoles) = 1;

[X,Y,T,AUC_MNE,OPTROCPT] = perfcurve(labels,Q_MNE_hat,1) ;
figure()
plot( X , Y)
title("  MNE")
hold on
scatter(OPTROCPT(1),OPTROCPT(2))
T = T(X==OPTROCPT(1) & Y==OPTROCPT(2));
predicted_MNE = (Q_MNE_hat>= T);
AUC_MNE
err_MNE = sqrt( sum( (labels- double(predicted_MNE)').^2,'all'))/length(Q_MNE)/3
acc_MNE = 1- err_MNE
%% WMNE
labels = zeros(1,length(Q_MNE)/3);
labels(index_dipoles) = 1;

[X,Y,T,AUC_WMNE,OPTROCPT] = perfcurve(labels,Q_WMNE_hat,1) ;
figure()
plot( X , Y)
title(" WMNE")
hold on
scatter(OPTROCPT(1),OPTROCPT(2))
T = T(X==OPTROCPT(1) & Y==OPTROCPT(2));
predicted_WMNE = (Q_WMNE_hat>= T);
AUC_WMNE
err_WMNE = sqrt( sum( (labels- double(predicted_WMNE)').^2,'all'))/1317
acc_WMNE = 1- err_WMNE
%% loreta 
d=1;

A1=zeros(length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    for j=1:length(Q_MNE)/3
        a=sqrt(sum((LocMat(:,i)-LocMat(:,j)).^2,'all'));
        if a==d
            A1(i,j)=1/6;
        end
    end
end
A0=inv(diag(A1*ones(length(Q_MNE)/3,1)))*A1;
A=kron(A0,eye(3));
B=(6/d^2)*(A-eye(length(Q_MNE)));

G=GainMat;
w=zeros(length(Q_MNE)/3,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
   w(i,i)=sqrt(sum( G(:,(i-1)*3+1:i*3)*G(:,(i-1)*3+1:i*3).','all'));
end
W=kron(w,eye(3));
W_loreta=W*B.'*B*W;

Q_Loreta=inv(W_loreta.'*W_loreta)*G.'*inv(G*(inv(W_loreta.'*W_loreta))*G.'+alpha*eye(21))*mean_peak;
Q_Loreta_hat=zeros(1,length(Q_MNE)/3);
for i=1:length(Q_MNE)/3
    Q_Loreta_hat(1,i)=sqrt(sum(Q_Loreta((i-1)*3+1:i*3).^2,'all'));
end
[Q index]=sort( Q_Loreta_hat ,'ascend');
labels = zeros(1,1317);
labels(index_dipoles) = 1;
[X,Y,T,AUC_Loreta,OPTROCPT] = perfcurve(labels,Q_Loreta_hat,1) ;
figure()
plot( X , Y)
title(" Loreta ")
hold on
scatter(OPTROCPT(1),OPTROCPT(2))
T = T(X==OPTROCPT(1) & Y==OPTROCPT(2));
predicted_Loreta = (Q_Loreta_hat>= T);

err_Loreta = sqrt( sum( (labels- double(predicted_Loreta)').^2,'all'))/1317
acc_Loreta = 1- err_Loreta