clc;
close all;
clear all;
rng(4)
S=10;
safe_distance = 0.2;
StartAgents = rand(S,2)*8 - 4;
% safe distance checking
distances = pdist(StartAgents); 
while any(distances < safe_distance)    
    StartAgents = rand(S, 2) * 8 - 4; 
    distances = pdist(StartAgents); % Regenerate positions
end
MovingAgents = StartAgents;

u_b = zeros(S,2);

% left_boundary = MovingAgents + [3.5,0];
% left_boundary = left_boundary(:,1) <= 0;
% u_b(:,1) = u_b(:,1) + 0.02*left_boundary;
% right_boundary = MovingAgents + [-3.5,0];
% right_boundary = right_boundary(:,1) >= 0;
% u_b(:,1) = u_b(:,1) - 0.02*right_boundary;
% upper_boundary = MovingAgents + [0,-3.5];
% upper_boundary = upper_boundary(:,2) >= 0;
% u_b(:,2) = u_b(:,2) - 0.02*upper_boundary;
% lower_boundary = MovingAgents + [0,3.5];
% lower_boundary = lower_boundary(:,2) <= 0;
% u_b(:,2) = u_b(:,2) + 0.02*lower_boundary;

boundary_offset = 3.5;

boundary_check = @(coord, dim, op) op*(coord(:,dim) + op*boundary_offset) <= 0;

u_b(:,1) = u_b(:,1) + 0.02*boundary_check(MovingAgents, 1, 1) - 0.02*boundary_check(MovingAgents, 1, -1);
u_b(:,2) = u_b(:,2) - 0.02*boundary_check(MovingAgents, 2, -1) + 0.02*boundary_check(MovingAgents, 2, 1);