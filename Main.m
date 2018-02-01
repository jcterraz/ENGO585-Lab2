% Author: Juan Carlos Terrazas Borbon
% Last Update: 2018-01-31
% Course: ENGO 585
% Lab: 2

% ---------------------Purpose of Code-------------------------------------
% The purpose of this code is to able to perform the first task required to
% perfom in the lab handout which involves Parametric Least Squares adjustment 

% Clear variables, close figure and command
clc
clear all
close all

% Read the file with data and extract its data
ranges = load('Lab2data.txt');

% Given Target Coordinates in meters
 targets = [0,0; 100,0; 100,100; 0,100];

%% Task 1: Batch Parametric Least Squares
est_coords = [50, 50];
P = diag(ones(4,1));

% 1.a Compute 2-D Solution for each epoch----------------------------------
x_hat_1_a = zeros(150,2);
for i = 1 : length(ranges)
    thres = 0;
    while thres == 0
        % Obtain the A matrix
        A = zeros(4,2); 
        for j = 1 : 4
            A(j, 1) = (est_coords(1) - targets(j, 1)) / ranges(i, j + 1);
            A(j, 2) = (est_coords(2) - targets(j, 2)) / ranges(i, j + 1);
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
        if abs(delta(1)) < 0.0001 && abs(delta(2)) < 0.0001
           thres = 1;     
        else
           est_coords = [est_coords(1) + delta(1),est_coords(2) + delta(2)];
        end
    end
    x_hat_1_a(i,:) = [est_coords(1) + delta(1),est_coords(2) + delta(2)];
    est_coords = [est_coords(1) + delta(1),est_coords(2) + delta(2)];
end

figure
plot(targets(:,1), targets(:,2),'o')
hold on
plot(x_hat_1_a(:,1), x_hat_1_a(:,2),'*');
hold off
title('Task 1.a 2-D Solution for each epoch')
xlabel('X Coordinates (meters)')
ylabel('Y Coordinates (meters)')
legend('Targets Points', 'Epoch Solution')

% 1.b Batch Solution-------------------------------------------------------
est_coords = [50, 50];
A_b = zeros(50*4,2);
w_b= zeros(50*4,1);
P_b = diag(ones(50*4,1));
thres = 0;
while thres == 0
    for i = 1: 50
        % Obtain the A matrix
        for j = 1 : 4
            A_b((i-1)*4+j, 1) = (est_coords(1) - targets(j, 1)) / ranges(i, j + 1);
            A_b((i-1)*4+j, 2) = (est_coords(2) - targets(j, 2)) / ranges(i, j + 1);
        end

        % Compute w Matrix
        for j = 1 : 4 
            w_b((i-1)*4+j, 1) = sqrt((targets(j, 1) - est_coords(1))^2 + ...
                (targets(j, 2) - est_coords(2))^2) - ranges(i, j + 1);
        end
    end
    % Compute N Matrix and obtain the delta values
    N = A_b' * P_b * A_b;
    delta = -1 * inv(N) * A_b' * P_b * w_b;
    
    % Check for delta and if threshold passes obtain coordinates
    if abs(delta(1)) < 0.0001 && abs(delta(2)) < 0.0001
        thres = 1;
        x_hat_1_b = [est_coords(1) + delta(1),est_coords(2) + delta(2)];
    else
        est_coords = [est_coords(1) + delta(1),est_coords(2) + delta(2)];
    end
end

% 1.c Plot Residuals-------------------------------------------------------
res = A * delta + w;

figure
plot(1:1:50, res(1:4:200,1))
title('Residuals from Target 1 (0,0)')
xlabel('Number of Measurement')
ylabel('Residual (meters)')

Apos =(res' * P * res)/(50-4);%aposteriori

C_x = inv(A' * P * A);%covariance matrix

% 1.d Error Ellipse--------------------------------------------------------
a = sqrt(0.5 * (C_x(1,1) + C_x(2,2)) + sqrt((1/4) * (C_x(1,1) - C_x(2,2))^2 ...
    + C_x(1,2)^2))* 2.45;

b = sqrt(0.5 * (C_x(1,1) + C_x(2,2)) - sqrt((1/4) * (C_x(1,1) - C_x(2,2))^2 ...
    + C_x(1,2)^2)) * 2.45;

azimuth = (1/2)* atan((2 * C_x(1,2))/(C_x(1,1) - C_x(2,2)));

t = linspace(0, 2*pi);
Theta = (90-azimuth) * pi/180;

figure
plot(x_hat_1_b(1),x_hat_1_b(2), 'o')
hold on
plot(x_hat_1_a(1:50,1),x_hat_1_a(1:50,2), '*')
plot(x_hat_1_b(1) + a*cos(Theta*t), x_hat_1_b(2) + b*sin(Theta*t))
legend('Batch Solution', 'Each Point Solution', 'Ellipse')
xlabel('East (meters)')
ylabel('North (meters)')

%% Task 2: Summation of Normals and Sequential LS
% 2.a Summation of Normal of 1.b


% 2.b Sequential Least Squares of 1.b


%% Task 3: Kalman Filtering
% 3.a Sequential Least Squares of whole data

% 3.b Kalman Filter of whole data


