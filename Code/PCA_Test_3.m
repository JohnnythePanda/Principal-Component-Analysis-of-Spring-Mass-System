close all; clear all; clc;

%%

load('cam1_3.mat');

close all;
vidFrames1_3 = vidFrames1_3(:,:,:,8:237);

x_pos = 330;
y_pos = 290;

numFrames = size(vidFrames1_3,4);
for j = 1:numFrames
    X = rgb2gray(vidFrames1_3(:,:,:,j));
    
    X(:, 1:x_pos-10) = 0;
    X(:, x_pos+20:640) = 0;
    
    X(1:y_pos-50, :) = 0;
    X(y_pos+10:480, :) = 0;
    
    [Max, Index] = max(X(:));
    
    [y_pos, x_pos] = ind2sub(size(X), Index);
    x_1(j) = x_pos;
    y_1(j) = y_pos;
%     imshow(X);
%     hold on 
%     plot(x_1(j), y_1(j), 'o', 'MarkerSize', 10);
%     drawnow
% %     

end

%% 

load('cam2_3.mat');
close all;
vidFrames2_3 = vidFrames2_3(:,:,:,36:265);

numFrames2 = size(vidFrames2_3,4);

x_pos = 360;
y_pos = 280;

for j = 1:numFrames2
    X2 = rgb2gray(vidFrames2_3(:,:,:,j));
    
     X2(:, 1:x_pos-10) = 0;
     X2(:, x_pos+10:640) = 0;
     
     X2(1:y_pos-10,:) = 0;
     X2(y_pos+10:480, :) = 0;
    
    [Max, Index] = max(X2(:));
    
    [y_pos, x_pos] = ind2sub(size(X2), Index);
    x_2(j) = x_pos;
    y_2(j) = y_pos;
%     imshow(X2);
%     hold on 
%     plot(x_2(j), y_2(j), 'o', 'MarkerSize', 10);
%     drawnow
%     

end

%%

load('cam3_3.mat');
close all;
vidFrames3_3 = vidFrames3_3(:,:,:,1:end-7);

numFrames3 = size(vidFrames3_3,4);

x_pos = 380;
y_pos = 250;

for j = 1:numFrames3
    X3 = rgb2gray(vidFrames3_3(:,:,:,j));
    
    X3(1:y_pos-20, :) = 0;
    X3(y_pos+20:480, :) = 0;
    
    X3(:, 1:x_pos-20) = 0;
    X3(:, x_pos+20:640) = 0;
    
    [Max, Index] = max(X3(:));
    
    [y_pos, x_pos] = ind2sub(size(X3), Index);
    x_3(j) = x_pos;
    y_3(j) = y_pos;
%     imshow(X3);
%     hold on 
%     plot(x_3(j), y_3(j), 'o', 'MarkerSize', 10);
%     drawnow
    

end

%% PCA 

%% PCA 

X = [x_1; y_1; x_2; y_2; x_3; y_3]; %% Stacks all row vectors into single matrix

figure(2);
plot(1:numFrames, X); %% Plots raw data, the coordinates of the mass from each camera.
set(gca,'Fontsize',14)
title("Test 3: Original Data")
xlabel("Time(Frames)")
ylabel("Position");
saveas(gcf,'Test 3 Original Data.jpg')


[m,n] = size(X); %% Stores the number of rows as m and columns as n.

mn = mean(X,2); %% Creates a 6x1 matrix containing the mean of each row of X
X = X - repmat(mn, 1, n); %% Subtracts the mean of each row from X. 
A = X/sqrt(n-1); %% Rescales X so we can use SVD on it. 

[U,S,V] = svd(A, 'econ'); %% Performs SVD on A. 

figure(3);
plot(V(:, 1:2)); %% Plots the time-series data
title("Test 3: Time-Series Data")
set(gca,'Fontsize',14)
xlabel("Time(Frames)")
ylabel("Position");

figure(4);
plot((diag(S).^2/sum(diag(S).^2)*100), 'o', 'Markersize', 10); %% Plots percentage
set(gca,'Fontsize',14)                                         %values of each sigma
title("Test 3: Energies of Principal Components")
xlabel("Principal Component")
ylabel("Energy(%)");
saveas(gcf,'Test 3 Energies of Singular Values.jpg')

U_star = U(:,1:2); %% A 6x2 matrix, each column is a principal component. 
projection = U_star'*(X); %% Projection of data using 2 principal components

figure(5);
plot(1:n,projection); %% Plots the projection data using 2 principal components. 
set(gca,'Fontsize',14)
title("Test 3: Projected Data")
xlabel("Time(Frames)")
ylabel("Amplitude");
legend('PCA Mode 1', 'PCA Mode 2', 'Location', 'SouthEast')
saveas(gcf,'Test 3 Projected Data.jpg')

