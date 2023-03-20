% Title: Global Peak Finding based on Distributed Gaussian Process and Active Sensing
% Author: Weijie Qi, Yunru Qu
% Date  : 01-2023

%%
clc;
close all;
clear all;
InputSpace = {linspace(-4, 4, 100); linspace(-4, 4, 100)};

%% The parameters here can be modified by the users

S = 7; % The number of stationary sensors is S*S
A = 10; % The number of moving agents is A
r = 1.8; % The range of communication with surrounding sensors

%% Show the environment

ShowEnvironment3D(InputSpace);
ShowTopDownView(InputSpace);

%% Compute the eigenfunctions and eigenvalues of Inputspace
[Eigenfunctions, Eigenvalues] = ComputeEigen(InputSpace);

%% Stationary sensors

[temp_x, temp_y] = meshgrid(linspace(-3.5, 3.5, S));
StationarySensors = zeros(S*S, 2);
StationarySensors(:, 1) = temp_x(:);
StationarySensors(:, 2) = temp_y(:);
Adj = BuildAdj(StationarySensors, r); % Adjacency matrix
y_s = f(temp_x, temp_y);
y_s = y_s(:);
% for i = 1:S*S
%     @TODO
% end

%% Moving agents

% Spawn the agents at random initial locations
AgentLocation = rand(A,2)*8 - 4;

% for d = 1:2
