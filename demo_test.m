clc;
clear; 
%% load data
img = imread('testimage\butterfly.jpg');
T= im2double(img);
Nway =size(T);

%% set missing indexes
ObsRatio=0.1;
Omega = randperm(prod(Nway));
Omega = Omega(1:round((ObsRatio)*prod(Nway)));
O = zeros(Nway);
O(Omega) = 1;
y=T.*O;
known=find(y);

%% set optimization paremeters
rho=0.0005; %for anisotropio TV 
w=[4,4,0];
[X_as] = Smoothlowrank_TVas( y, known, rho,y,w);

 rho=0.002; %for isotropio TV 
  w=[2,2,0];
[X_is] = Smoothlowrank_TVis( y, known, rho,y,w);

figure(1),
subplot(2,2,1),imshow(T),title('orginal image');
subplot(2,2,2),imshow(y),title('missiong image');
subplot(2,2,3),imshow(X_as),title('X as');
subplot(2,2,4),imshow(X_is),title('X is');

