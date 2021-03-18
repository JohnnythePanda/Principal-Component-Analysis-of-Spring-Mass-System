close all; clear all; clc;

%%

load('cam1_1.mat'); %% Load in video in width by height by rgb by frame matrix

x_pos = 319; %% Initial x location of flashlight/area I am tracking
y_pos = 225; %% Initial y location of flashlight/area I am tracking

vidFrames1_1 = vidFrames1_1(:,:,:,7:223); %% Trims the video to contain only
                                          % frames 7 to 223 of the video. 
numFrames = size(vidFrames1_1,4); %% Number of frames in video. 
for j = 1:numFrames
    X = rgb2gray(vidFrames1_1(:,:,:,j)); %% Converts the j-th frame to grayscale
    
    X(:, 1:x_pos-20) = 0; %% Makes the frame zero outside 20 pixels region
    X(:, x_pos+20:640) = 0; %% from the location of the tracked area. 
    
    X(1:y_pos-20, :) = 0;
    X(y_pos+20:480, :) = 0;
    
    [Max, Index] = max(X(:)); %% Finds the pixel with highest value
    
    [y_pos, x_pos] = ind2sub(size(X), Index); %% Finds the indexes of the pixel
                                              %with highest value
                                              
    x_1(j) = x_pos; %% Stores the x-coordinate of mass in j-th frame in vector           
    y_1(j) = y_pos; %% Stores the x-coordinate of mass in j-th frame in vector
    
%     imshow(X);  %% Plots the filtered grayscale j-th frame of the video.  
%     hold on 
%     plot(x_1(j), y_1(j), 'o', 'MarkerSize', 10); %% Plots the coordinates being stored
%                                                  % on the j-th frame. 
%     drawnow
    

end

%% 

load('cam2_1.mat');

vidFrames2_1 = vidFrames2_1(:,:,:,17:233);

numFrames2 = size(vidFrames2_1,4);

x_pos = 278;
y_pos = 122;

for j = 1:numFrames2
    X2 = rgb2gray(vidFrames2_1(:,:,:,j));
    
     X2(:, 1:x_pos-20) = 0;
     X2(:, x_pos+20:640) = 0;
     
     X2(1:y_pos-30,:) = 0;
     X2(y_pos+30:480, :) = 0;
    
    [Max, Index] = max(X2(:));
    
    [y_pos, x_pos] = ind2sub(size(X2), Index);
    x_2(j) = x_pos;
    y_2(j) = y_pos;
    
%     imshow(X2);
%     hold on 
%     plot(x_2(j), y_2(j), 'o', 'MarkerSize', 10);
%     drawnow
    

end

%%

load('cam3_1.mat');

vidFrames3_1 = vidFrames3_1(:,:,:,6:222);

numFrames2 = size(vidFrames3_1,4);

x_pos = 295;
y_pos = 267;

for j = 1:numFrames2
    X3 = rgb2gray(vidFrames3_1(:,:,:,j));
    
    X3(1:y_pos-20, :) = 0;
    X3(y_pos+20:480, :) = 0;
    
    X3(:, 1:x_pos-30) = 0;
    X3(:, x_pos+30:640) = 0;
    
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

X = [x_1; y_1; x_2; y_2; x_3; y_3]; %% Stacks all row vectors into single matrix

figure(2);
plot(1:numFrames, X); %% Plots raw data, the coordinates of the mass from each camera.
set(gca,'Fontsize',14)
title("Test 1: Original Data")
xlabel("Time(Frames)")
ylabel("Position");
saveas(gcf,'Test 1 Original Data.jpg')


[m,n] = size(X); %% Stores the number of rows as m and columns as n.

mn = mean(X,2); %% Creates a 6x1 matrix containing the mean of each row of X
X = X - repmat(mn, 1, n); %% Subtracts the mean of each row from X. 
A = X/sqrt(n-1); %% Rescales X so we can use SVD on it. 

[U,S,V] = svd(A, 'econ'); %% Performs SVD on A. 

figure(3);
plot(V(:, 1:2)); %% Plots the time-series data
title("Test 1: Time-Series Data")
set(gca,'Fontsize',14)
xlabel("Time(Frames)")  
ylabel("Position");

figure(4);
plot((diag(S).^2/sum(diag(S).^2)*100), 'o', 'Markersize', 10); %% Plots percentage
set(gca,'Fontsize',14)                                         %values of each sigma
title("Test 1: Energies of Principal Components")
xlabel("Principal Component")
ylabel("Energy(%)");
saveas(gcf,'Test 1 Energies of Singular Values.jpg')

U_star = U(:,1:2); %% A 6x2 matrix, each column is a principal component. 
projection = U_star'*(X); %% Projection of data using 2 principal components

figure(5);
plot(1:n,projection); %% Plots the projection data using 2 principal components. 
set(gca,'Fontsize',14)
title("Test 1: Projected Data")
xlabel("Time(Frames)")
ylabel("Amplitude");
legend('PCA Mode 1', 'PCA Mode 2', 'Location', 'SouthEast')
saveas(gcf,'Test 1 Projected Data.jpg')

