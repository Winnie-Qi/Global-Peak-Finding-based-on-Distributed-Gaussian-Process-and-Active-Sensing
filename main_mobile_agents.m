% Title: Global Peak Finding based on Distributed Gaussian Process and Active Sensing
% Author: Weijie Qi, Yunru Qu
% Date  : 07-2023

clc;
close all;
clear all;

% The parameters here can be modified by the users
S = 6; % The number of moving agents is S
r = 1.0; % The range of communication with nearby agents
a = 100; % test points density
E = 100;

% hyperparameters
sigma = 0.0001;
k_max = 720; % max iteration times
SHOW = [1,6]; % show the result of the Sth sensor, please note that the more agents you select, the slower the results appear
SensingPeriod = 30; % The number of steps between obtaining new measurements
gamma = 0.3;
v = 0.01;
l = 0.01;
% rng('default'); % If you want to generate new initial positions every time, comment out this line
rng(2)

figure(1)
InputSpace_test = {linspace(-4, 4, a); linspace(-4, 4, a)};

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
MovingAgents = StartAgents;

% hold on
% h_drones = plot3(MovingAgents(:,1), MovingAgents(:,2), z, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
% for i = 1:S
%     plot3([MovingAgents(i,1), MovingAgents(i,1)], [MovingAgents(i,2), MovingAgents(i,2)], [0, z(i)], '--k');
% end
% [~,~] = BuildAdj(MovingAgents, r);

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

% f_E = cell(S,1);
% Pi_E = cell(S,1);

ExplorationFlag = ones(1,S);
L = 20*ones(1,S);
reachGoal = zeros(1,S);
steps = inf(1,S);
k = 0; m = 1;
boundary_offset = 2.9;
boundary_check = @(coord, dim, op) op*(coord(:,dim) + op*boundary_offset) <= 0;

v = zeros(S,2);

% main loop
while k <= k_max        
    [distances,tran] = BuildAdj(MovingAgents, r, InputSpace_test, S);    
    u_b = zeros(S,2); % boundary condition
    u_c = zeros(S,2); % collision condition
%     u_g = zeros(S,2); % goal oriented
    
    % boundary check
    u_b(:,1) = u_b(:,1) + 0.02*boundary_check(MovingAgents, 1, 1) - 0.02*boundary_check(MovingAgents, 1, -1);
    u_b(:,2) = u_b(:,2) - 0.02*boundary_check(MovingAgents, 2, -1) + 0.02*boundary_check(MovingAgents, 2, 1);

    % collision check
    for i=1:S
        if ExplorationFlag
            for j=1:S
                u_c(i,:) = u_c(i,:) + distances(i,j) * 0.003./(MovingAgents(i,:)-MovingAgents(j,:)+0.00001); 
                u_c = min(max(u_c, -0.05), 0.05);
            end
        end
    end
        
    if mod(k,SensingPeriod) == 0        

        % Measurement update
        y_s = f(MovingAgents(:,1), MovingAgents(:,2)) + sigma*randn(size(MovingAgents,1),1); % get new measurements
        v = zeros(S,2);
%         [distances,tran] = BuildAdj(MovingAgents, r); % get new adjacency matrix % 这个函数里面注释掉了一些       
        
        for n = 1:S     
            % update alpha and beta
            clear phi tmp Pi_E
            phi = Find_Eigenfunctions_by_LUT(MovingAgents(n,:),E);          
            tmp = phi' * phi;
            alpha(n,:) = (m-1)/m * alpha(n,:) + 1/m * tmp(:)';
            beta(n,:) = (m-1)/m * beta(n,:) + 1/m * (phi' * y_s(n))'; 
%             alpha(n,:) = alpha(n,:) + tmp(:)';
%             beta(n,:) = beta(n,:) + (phi' * y_s(n))';             
            
            % compute acceleration or velocity
            if ExplorationFlag(n)                
%                 [goto_x,goto_y] = findMaxEdge(MovingAgents(n,:),reshape(diag(Pi_E{n}),a,a),1);
                Pi_E = K - PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * reshape(alpha(n,:),E,E)*LAMBDA*PHI';
                [goto_x,goto_y,ExplorationFlag(n),~] = findMaxEdge(MovingAgents(n,:),reshape(diag(Pi_E),a,a),1);
                v(n,:) = 0.03 * [goto_x,goto_y];
                 
            else
%                 f_E{n} = PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(n,:)';
%                 [goto_x,goto_y] = findMaxEdge(MovingAgents(n,:),reshape(f_E{n},a,a),0);
                f_E = PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(n,:)';
                [goto_x,goto_y,~,step] = findMaxEdge(MovingAgents(n,:),reshape(f_E,a,a),0);
                v(n,:) = 0.03 * [goto_x,goto_y];
%                 disp(n)
%                 disp(v(n,:))
                if step
                    steps(n) = step;
                end                
            end
            
        end        
        
        m = m + 1;       
    end
    for n = 1:S
        if ExplorationFlag
%             f_E{n} = PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(n,:)';
        end
%         Pi_E{n} = K - PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * reshape(alpha(n,:),E,E)*LAMBDA*PHI';         
    end    
            
    % 临时
%     f_E{SHOW} = PHI * pinv(reshape(alpha(SHOW,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(SHOW,:)';
%     Pi_E{SHOW} = K - PHI * pinv(reshape(alpha(SHOW,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * reshape(alpha(SHOW,:),E,E)*LAMBDA*PHI';

    if mod(k,SensingPeriod) == 0        
        for n = 1:S
            f_E = PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(n,:)';                       
            subplot(4,S,2*S+n);
            pcolor(tmp_xx,tmp_yy,reshape(f_E,a,a));            
            clear f_E            
            shading flat
            Pi_E = K - PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * reshape(alpha(n,:),E,E)*LAMBDA*PHI'; 
            subplot(4,S,3*S+n);
            pcolor(tmp_xx,tmp_yy,reshape(diag(Pi_E),a,a));            
            shading flat            
            clear Pi_E
            pause(0.05);
        end

    end
    
    % update alpha and beta through communication
    alpha = alpha - gamma*tran*alpha;
    beta = beta - gamma*tran*beta;    
    
    % moving
    v = v + u_b + u_c;
%     v = min(max(v, -0.05), 0.05);
%     q = 2;
% %     disp(k);
%     disp('u_b')
%     disp(u_b(q,:))
%     disp('u_c')
%     disp(u_c(q,:))
%     disp('v')
%     disp(v(q,:));
    MovingAgents = MovingAgents + (~reachGoal' .* v);
%     MovingAgents = MovingAgents + v;
    steps = steps - 1;
    reachGoal = reachGoal + (steps == 0);
    
    disp(ExplorationFlag)
%     disp(reachGoal)
%     disp(steps)
%     disp('MovingAgents')
%     disp(MovingAgents(q,:))
%     MovingAgents(:,1) = MovingAgents(:,1)+0.1; % moving
    % 在这里可以画一次激光效果图
    
%     delete(h_drones);        
%     h_drones = plot3(MovingAgents(:,1), MovingAgents(:,2), z, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
    k = k + 1;
    disp(k)
    pause(0.01);    
end
