% Title: Global Peak Finding based on Distributed Gaussian Process and Active Sensing
% Author: Weijie Qi, Yunru Qu
% Date  : 07-2023

clc;
close all;
clear all;

% The parameters here can be modified by the users
S = 6; % The number of moving agents is S
rng(1) % If you want to generate new initial positions every time, comment out this line

% hyperparameters
r = 1.0; % The range of communication with nearby agents
a = 100; % test points density
E = 120;
sigma = 0.0001;
k_max = 500; % max iteration times
% SHOW = [1,6]; % show the result of the Sth sensor, please note that the more agents you select, the slower the results appear
SensingPeriod = 30; % The number of steps between obtaining new measurements
% gamma = 0.3;
gamma = 0.00001;
% gamma = 0.01;
v = 0.01;
l = 0.01;

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

% Intialize alpha and beta
alpha = zeros(S,E*E); % (1,E*E) every sensor
beta = zeros(S,E); % (1,E) every sensor


ExplorationFlag = ones(1,S);
L = 20*ones(1,S);
reachGoal = zeros(1,S);
update = ones(1,S);
k = 0; m = 1;
boundary_offset = 3;
boundary_check = @(coord, dim, op) op*(coord(:,dim) + op*boundary_offset) <= 0;

v = zeros(S,2);

% main loop
while k <= k_max        
    ray = mod(k,SensingPeriod);
    [distances,tran] = BuildAdj(MovingAgents, r, InputSpace_test, S, ray);    
    u_b = zeros(S,2); % boundary condition
    u_c = zeros(S,2); % collision condition
    
    % boundary check
    u_b(:,1) = u_b(:,1) + 0.01*boundary_check(MovingAgents, 1, 1) - 0.01*boundary_check(MovingAgents, 1, -1);
    u_b(:,2) = u_b(:,2) - 0.01*boundary_check(MovingAgents, 2, -1) + 0.01*boundary_check(MovingAgents, 2, 1);

    % collision check
    for i=1:S        
        for j=1:S
            u_c(i,:) = u_c(i,:) + distances(i,j) * 0.003./(MovingAgents(i,:)-MovingAgents(j,:)+0.00001);
            if ExplorationFlag(i)
                u_c(i,:) = min(max(u_c(i,:), -0.1), 0.1);
            else
                u_c(i,:) = min(max(u_c(i,:), -0.003), 0.003);
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
            if update(n) % It cannot be updated too densely because the accuracy of LUT is not high enough
                alpha(n,:) = (m-1)/m * alpha(n,:) + 1/m * tmp(:)';
                beta(n,:) = (m-1)/m * beta(n,:) + 1/m * (phi' * y_s(n))'; 
%             else
%                 beta(n,:) = (m-1)/m * beta(n,:) + 1/m * (phi' * y_s(n))';
%             elseif ~reachGoal(n)
%                 alpha(n,:) = (m + 5*floor((k - 570)/30)-1)/(m + 5*floor((k - 570)/30)) * alpha(n,:) + 1/(m + 5*floor((k - 570)/30)) * tmp(:)';
%                 beta(n,:) = (m + 5*floor((k - 570)/30)-1)/(m + 5*floor((k - 570)/30)) * beta(n,:) + 1/(m + 5*floor((k - 570)/30)) * (phi' * y_s(n))';
            end
%             alpha(n,:) = alpha(n,:) + tmp(:)';
%             beta(n,:) = beta(n,:) + (phi' * y_s(n))';
%             
            % compute acceleration or velocity
            if ExplorationFlag(n)                
                Pi_E = K - PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * reshape(alpha(n,:),E,E)*LAMBDA*PHI';
                [goto_x,goto_y,ExplorationFlag(n),~] = findMaxEdge(MovingAgents(n,:),reshape(diag(Pi_E),a,a),1);
                v(n,:) = 0.2 * [goto_x,goto_y]/SensingPeriod;
                 
            else
                f_E = PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(n,:)';
                [goto_x,goto_y,~,reach] = findMaxEdge(MovingAgents(n,:),reshape(f_E,a,a),0);
                v(n,:) = 0.16 * [goto_x,goto_y]/SensingPeriod;
                if reach
                    reachGoal(n) = reach;
                end                
            end
            
        end        
        
        m = m + 1;       
    end

 % plot
 labels = cellstr(num2str((1:S)', 'agent %d'));
    if mod(k,SensingPeriod) == 0        
        for n = 1:S            
            
            f_E = PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(n,:)';                       
            subplot(4,S,2*S+n);
            pcolor(tmp_xx,tmp_yy,reshape(f_E,a,a)');            
            clear f_E            
            shading flat            
            Pi_E = K - PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * reshape(alpha(n,:),E,E)*LAMBDA*PHI'; 
            Pi_E = diag(Pi_E);
            if sum(Pi_E(:)) < 70000
                update(n) = 0;
            end
            subplot(4,S,3*S+n);
            pcolor(tmp_xx,tmp_yy,reshape(Pi_E,a,a));            
            shading flat            
            clear Pi_E
            
            % text
            subplot(4, S, 2*S+1);
            text(-7, 0.5, 'estimate', 'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');            
            subplot(4, S, 3*S+1);
            text(-7, 0.5, 'uncertainty', 'FontSize', 7, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');            
            subplot(4, S, 2*S+n); 
            text(0.5, 3.7, labels{n}, 'FontSize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            subplot(4, S, 3*S+n)
            if ExplorationFlag(n)
                text(0.5, -7, 'explorating...', 'Color', 'red', 'FontSize', 7, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            elseif ~reachGoal(n)
                text(0.5, -7, 'goal targeting...', 'Color', 'blue', 'FontSize', 7, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            else
                text(0.5, -7, 'reached goal', 'Color', 'black', 'FontSize', 7, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
            
            
            pause(0.05);
        end

    end
    
    % update alpha and beta through communication
    tran(~update, :) = 0;
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
%     MovingAgents = MovingAgents + (~reachGoal' .* v);
%     MovingAgents = MovingAgents + (reachGoal>0)'.* v * 0.25;
%     MovingAgents = MovingAgents + ~reachGoal' .* v;
%     MovingAgents = MovingAgents + (reachGoal' .* v)*1/2;
    MovingAgents = MovingAgents + (ExplorationFlag' .* v);
    MovingAgents = MovingAgents + ((~ExplorationFlag& ~reachGoal)' .* v)*1/2;
    MovingAgents = MovingAgents + (reachGoal' .* v)*1/3;
%     steps = steps - 1;
%     reachGoal = reachGoal + (steps == 0);
    
%     disp(ExplorationFlag)
    disp(update)
%     disp(steps)
%     disp('MovingAgents')
%     disp(MovingAgents(q,:))
%     MovingAgents(:,1) = MovingAgents(:,1)+0.1; % moving
    
    k = k + 1;     
%     if k>500
%         gamma = 0;
%     else
    if k>100
        gamma = 0.1^(floor((k-100)/100) + 6);
    end

    disp(k)
    pause(0.01);    
end
