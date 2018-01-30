% Author: Juan Carlos Terrazas Borbon
% Last Update: 2018-01-29
% Course: ENGO 585
% Lab: 2

% ---------------------Purpose of Code-----------------------------------
% The purpose of this code is to able to perform the 3 task required to
% perfome in the lab handout which involves Least Squares adjustment and
% Kalman filter

% Clear variables, close figure and command
clc
clear all
close all

% Read the file with data and extract its data
ranges = load('Lab2data.txt');
% range.time = data(:,1);
% range.one = data(:,2); 
% range.two = data(:,3); 
% range.three = data(:,4); 
% range.four = data(:,5);
% 
% clear data

% Given Target Coordinates in meters
% target.one = [0,0];
% target.two = [100,0];
% target.three = [100,100];
% target.four = [0,100];
 targets = [0,0; 100,0; 100,100; 0,100];

%% Task 1: Batch Parametric Least Squares
est_coords = [0, 0];
P = diag(ones(4,1));

% 1.a Compute 2-D Solution for each epoch
x_hat = zeros(150,2);
for i = 1 : length(ranges)
    thres = 0;
    while thres == 0
        % Obtain the A matrix
        A = zeros(4,2); 
        for j = 1 : 4
            A(j, 1) = (targets(j, 1) - est_coords(1)) / ranges(i, j + 1);
            A(j, 2) = (targets(j, 2) - est_coords(2)) / ranges(i, j + 1);
        end
        
        % Compute w Matrix
        w= zeros(4,1);
        for j = 1 : 4 
            w(j, 1) = sqrt((targets(j, 1) - est_coords(1))^2 + ...
                (targets(j, 2) - est_coords(2))^2) - ranges(i, j + 1);
        end
        
        % Compute N Matrix and obtain the delta values
        N = A' * P * A;
        delta = -1 * inv(N) * A' * P * w;
        
        % check 
        if any(abs(delta) < 0.0001)
           thres = 1;     
        else
           est_coords = [est_coords(1) + delta(1),est_coords(2) + delta(2)];
        end
    end
    x_hat(i,:) = [est_coords(1) + delta(1),est_coords(2) + delta(2)];
end

%% Task 2: Summation of Normals and Sequential LS


%% Task 3: Kalman Filtering

