% Title: Global Peak Finding based on Distributed Gaussian Process and Active Sensing
% Author: Weijie Qi, Yunru Qu
% Date  : 07-2023

clc;
close all;
clear all;

% The parameters here can be modified by the users
S = 20; % The number of moving agents is S
r = 1.0; % The range of communication with nearby agents
a = 100; % test points density
E = 40;
sigma = 0.0001;
M = 1000; % iteration times
SHOW = 6; % show the result of the Sth sensor
gamma = 0.02;
rng('default'); % If you want to generate new initial positions every time, comment out this line

% Show the environment
InputSpace_test = {linspace(-4, 4, a); linspace(-4, 4, a)};
% ShowEnvironment3D(InputSpace);
% ShowTopDownView(InputSpace_test);
ShowEnvironment3D(InputSpace_test);

% Load pre-computed LUT 
load mat/InputSpace_(-4,4,100)_E_463.mat
global Eig_LUT InputSpace_LUT;
Eig_LUT = Eigenfunctions;
InputSpace_LUT = InputSpace;

% Find the eigenfunctions and eigenvalues of Inputspace_test
PHI = Find_Eigenfunctions_by_LUT(InputSpace_test,E);
LAMBDA = diag(Eigenvalues(1:E));

% Initialize the positions of the moving agents
safe_distance = 0.2;
MovingAgents = rand(S,2)*8 - 4;
% safe distance checking
distances = pdist(MovingAgents); 
while any(distances < safe_distance)    
    MovingAgents = rand(20, 2) * 8 - 4; 
    distances = pdist(MovingAgents); % Regenerate positions
end
z = ones(S,1) * 15;
% drone_img = imread('drone.jpg');
% for i = 1:S    
%     x_pos = MovingAgents(S,1);
%     y_pos = MovingAgents(S,2);
%     z_pos = z(i);
%     x_range = [x_pos-0.2, x_pos+0.2];
%     y_range = [y_pos-0.2, y_pos+0.2];
%     z_range = [z_pos, z_pos];
%     texture_data = repmat(drone_img, [1, 1, 3]);
%     surf(x_range, y_range, z_range, texture_data, 'FaceColor', 'texturemap', 'EdgeColor', 'none');
% end
h_drones = plot3(MovingAgents(:,1), MovingAgents(:,2), z, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
for t = 0:0.1:2
    delete(h_drones);
%     scatter3(MovingAgents(:,1), MovingAgents(:,2), z, 'filled', 'r');
    MovingAgents(:,1) = MovingAgents(:,1)-0.1;
    h_drones = plot3(MovingAgents(:,1), MovingAgents(:,2), z, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    pause(0.5);
end
% for i = 1:S
%     plot3([MovingAgents(i,1), MovingAgents(i,1)], [MovingAgents(i,2), MovingAgents(i,2)], [0, z(i)], '--k');
% end