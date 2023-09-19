% Title: Global Peak Finding based on Distributed Gaussian Process and Active Sensing
% Author: Weijie Qi, Yunru Qu
% Date  : 07-2023

clc;
close all;
clear all;

% The parameters here can be modified by the users
S = 10; % The number of moving agents is S
r = 1.0; % The range of communication with nearby agents
a = 100; % test points density
E = 100;

% hyperparameters
sigma = 0.0001;
k_max = 30; % max iteration times
SHOW = 8; % show the result of the Sth sensor
SensingPeriod = 5; % The number of steps between obtaining new measurements
gamma = 0.7;
v = 0.01;
l = 0.01;
% rng('default'); % If you want to generate new initial positions every time, comment out this line
rng(2)

% Show the environment
InputSpace_test = {linspace(-4, 4, a); linspace(-4, 4, a)};
% ShowEnvironment3D(InputSpace);
ShowTopDownView(InputSpace_test);
% ShowEnvironment3D(InputSpace_test);

% Load pre-computed LUT 
load mat/InputSpace_(-4,4,100)_E_463.mat
global Eig_LUT InputSpace_LUT;
Eig_LUT = Eigenfunctions;
InputSpace_LUT = InputSpace;

% Find the eigenfunctions and eigenvalues of Inputspace_test
[test_x, test_y] = meshgrid(InputSpace_test{1},InputSpace_test{2});
TestPoints = zeros(a*a, 2);
TestPoints(:, 1) = test_x(:)';
TestPoints(:, 2) = test_y(:)';
PHI = Find_Eigenfunctions_by_LUT(TestPoints,E);
LAMBDA = diag(Eigenvalues(1:E));

% kernel of test points
K = zeros(a*a);
[tmp_xx, tmp_yy] = meshgrid(linspace(-4, 4, a));
tmp_x = tmp_xx(:);
tmp_y = tmp_yy(:);
for i = 1:a*a
    for j = 1:a*a        
        distance = sqrt((tmp_x(i) - tmp_x(j))^2 + (tmp_y(i) - tmp_y(j))^2);
        K(i,j) = exp(-0.5*distance^2/v)*l+11;
    end
end

% Initialize the positions of the moving agents
safe_distance = 0.2;
StartAgents = rand(S,2)*8 - 4;
% safe distance checking
distances = pdist(StartAgents); 
while any(distances < safe_distance)    
    StartAgents = rand(S, 2) * 8 - 4; 
    distances = pdist(StartAgents); % Regenerate positions
end
z = ones(S,1) * 15;
% StartAgents = [[-1.7720,-0.6259];[-1.6920,-0.6259]];
% StartAgents = [[-1.7720,-0.6259];[-1.6920,-0.6259];[-1.6120,-0.6259];[-1.532,-0.6259]];
MovingAgents = StartAgents;

y_s = f(MovingAgents(:,1), MovingAgents(:,2)) + sigma*randn(size(MovingAgents,1),1);

if 0 % if you wannt to see the distributed result here
    G = Find_Eigenfunctions_by_LUT(MovingAgents,E);
    f_E_central = PHI * pinv(G'*G/S^2+sigma^2/S^2*pinv(LAMBDA)) * G'/S^2*y_s; 
    figure(3)
    pcolor(tmp_xx,tmp_yy,reshape(f_E_central,a,a));
    colorbar
    k = K - PHI * pinv(G'*G/S^2+sigma^2/S^2*pinv(LAMBDA)) * G'/S *G *LAMBDA*PHI'; 
    figure(4)
    pcolor(tmp_xx,tmp_yy,reshape(diag(k),a,a));
    colorbar
end

% figure(5)
% h_drones = plot3(StartAgents(:,1), StartAgents(:,2), z, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Intialize alpha and beta
alpha = zeros(S,E*E); % (1,E*E) every sensor
beta = zeros(S,E); % (1,E) every sensor

f_E = cell(S,1);
Pi_E = cell(S,1);

ExplorationFlag = 1;
k = 0; m = 1;
% main loop
while k <= k_max
    if mod(k,SensingPeriod) == 0        
        y_s = f(MovingAgents(:,1), MovingAgents(:,2)) + sigma*randn(size(MovingAgents,1),1); % get new measurements        
        [Adj,tran] = BuildAdj(MovingAgents, r); % get new adjacency matrix % 这个函数里面注释掉了一些        
        % update alpha and beta
        for n = 1:S     
            phi = Find_Eigenfunctions_by_LUT(MovingAgents(n,:),E);          
            tmp = phi' * phi;
%             alpha(n,:) = tmp(:);
%             beta(n,:) = phi' * y_s(n);
%             alpha(n,:) = (m-1)/m * alpha(n,:) + 1/m * tmp(:)';
%             beta(n,:) = (m-1)/m * beta(n,:) + 1/m * (phi' * y_s(n))';
            alpha(n,:) = alpha(n,:) + tmp(:)';
            beta(n,:) = beta(n,:) + (phi' * y_s(n))';
        end
        m = m + 1;
    end
    for n = 1:S
        if ExplorationFlag
            f_E{n} = PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(n,:)';
        end
%         Pi_E{n} = K - PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * reshape(alpha(n,:),E,E)*LAMBDA*PHI';         
    end    
            
    % 临时
%     f_E{SHOW} = PHI * pinv(reshape(alpha(SHOW,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(SHOW,:)';
%     Pi_E{SHOW} = K - PHI * pinv(reshape(alpha(SHOW,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * reshape(alpha(SHOW,:),E,E)*LAMBDA*PHI';

    if mod(k,1) == 0
%         f_E{SHOW} = PHI * pinv(reshape(alpha(SHOW,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(SHOW,:)';
        Pi_E{SHOW} = K - PHI * pinv(reshape(alpha(SHOW,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * reshape(alpha(SHOW,:),E,E)*LAMBDA*PHI';
%         figure(k+10)
%         pcolor(tmp_xx,tmp_yy,reshape(f_E{SHOW},a,a));
%         colorbar
%         shading flat
        figure(k+100)
        pcolor(tmp_xx,tmp_yy,reshape(diag(Pi_E{SHOW}),a,a));        
        colorbar
        shading flat
        hold on
        plot(MovingAgents(SHOW,1),MovingAgents(SHOW,2),'o');
    end
    
    % update alpha and beta through communication
    alpha = alpha - gamma*tran*alpha;
    beta = beta - gamma*tran*beta;    
    
    MovingAgents(:,1) = MovingAgents(:,1)+0.1; % moving
    % 在这里可以画一次激光效果图
    
%     delete(h_drones);        
%     h_drones = plot3(MovingAgents(:,1), MovingAgents(:,2), z, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    k = k + 1;    
    pause(0.1);    
end
% for i = 1:S
%     plot3([StartAgents(i,1), StartAgents(i,1)], [StartAgents(i,2), StartAgents(i,2)], [0, z(i)], '--k');
% end