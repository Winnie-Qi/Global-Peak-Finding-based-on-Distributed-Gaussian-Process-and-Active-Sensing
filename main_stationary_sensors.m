% Title: Global Peak Finding based on Distributed Gaussian Process and Active Sensing
% Author: Weijie Qi, Yunru Qu
% Date  : 01-2023

clc;
close all;
clear all;

%% The parameters here can be modified by the users

S = 6; % The number of stationary sensors is S*S£¬corresponds to traing points density
r = 1.8; % The range of communication with surrounding sensors
a = 50; % test points density
E = 50;
sigma = 0.1;
M = 100; % iteration times
gamma = 0.02;

%%

% Show the environment
InputSpace = {linspace(-4, 4, a); linspace(-4, 4, a)};
% ShowEnvironment3D(InputSpace);
% ShowTopDownView(InputSpace);

% Compute the eigenfunctions and eigenvalues of Inputspace
PHI = zeros(a^2,E);
i = 1;
for t = 1:a
    for s = 1:a
        for e = 1:E
            PHI(i,e) = 2*sin((e-1/2)*pi*InputSpace{1}(s)) * 2*sin((e-1/2)*pi*InputSpace{2}(t));
        end
        i = i+1;
    end
end

LAMBDA = zeros(E,E);
for e = 1:E
    LAMBDA(e,e) = 1/(e*pi-pi/2)^4;
end

% Initialize the stationary sensors
[temp_x, temp_y] = meshgrid(linspace(-3.5, 3.5, S));
StationarySensors = zeros(S*S, 2);
StationarySensors(:, 1) = temp_x(:);
StationarySensors(:, 2) = temp_y(:);
Adj = BuildAdj(StationarySensors, r); % Adjacency matrix
y_s = f(temp_x, temp_y) + sigma*randn(1,1);
y_s = y_s(:);

% Intialize alpha and beta
alpha = zeros(E,E);
beta = zeros(E,1);

while true
    f_E = PHI * inv(alpha+sigma^2/S^2*inv(LAMBDA)) * beta;
end