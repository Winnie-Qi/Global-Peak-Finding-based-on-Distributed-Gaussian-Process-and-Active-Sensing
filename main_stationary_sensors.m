% Title: Global Peak Finding based on Distributed Gaussian Process and Active Sensing
% Author: Weijie Qi, Yunru Qu
% Date  : 01-2023

clc;
close all;
clear all;

% The parameters here can be modified by the users

S = 10; % The number of stationary sensors is S*S，corresponds to traing points density
r = 1.0; % The range of communication with surrounding sensors
a = 100; % test points density
E = 40;
sigma = 0.0001;
M = 20; % iteration times
SHOW = 10; % show the result of the Sth sensor
gamma = 0.7;
v = 0.05;
l = 0.05;
interval = 1;

% Show the environment
InputSpace_test = {linspace(-4, 4, a); linspace(-4, 4, a)};
% ShowEnvironment3D(InputSpace);
ShowTopDownView(InputSpace_test);

% Load pre-computed LUT 
load mat/InputSpace_(-4,4,100)_E_463.mat
global Eig_LUT InputSpace_LUT;
Eig_LUT = Eigenfunctions;
InputSpace_LUT = InputSpace;

% Find the eigenfunctions and eigenvalues of Inputspace_test
PHI = Find_Eigenfunctions_by_LUT(InputSpace_test,E);
LAMBDA = diag(Eigenvalues(1:E));

% kernel of test points
K = zeros(a*a);
[tmp_x, tmp_y] = meshgrid(linspace(-4, 4, a));
tmp_x = tmp_x(:);
tmp_y = tmp_y(:);
for i = 1:a*a
    for j = 1:a*a        
        distance = sqrt((tmp_x(i) - tmp_x(j))^2 + (tmp_y(i) - tmp_y(j))^2);
        K(i,j) = exp(-0.5*distance^2/v)*l+6;
    end
end

% Initialize the stationary sensors
InputSpace_train = {linspace(-3.5, 3.5, S);linspace(-3.5, 3.5, S)};
[train_x, train_y] = meshgrid(linspace(-3.5, 3.5, S));
StationarySensors = zeros(S*S, 2);
train_x = train_x';train_y = train_y';
StationarySensors(:, 1) = train_x(:)';
StationarySensors(:, 2) = train_y(:)';
Adj = BuildAdj(StationarySensors, r); % Adjacency matrix
y_s = f(train_x, train_y) + sigma*randn(S,S);
y_s = y_s(:);

if 1 % if you wannt to see the distributed result here  
    G = Find_Eigenfunctions_by_LUT(InputSpace_train,E);
    f_E_central = PHI * pinv(G'*G/S^2+sigma^2/S^2*pinv(LAMBDA)) * G'/S^2*y_s; 
    figure(3)
    pcolor(reshape(f_E_central,a,a));
    colorbar
end

% Intialize alpha and beta
alpha = zeros(S*S,E*E); % 每个sensor是(1,E*E)
beta = zeros(S*S,E); % 每个sensor是(1,E)
for i = 1:S
    for j = 1:S
        n = (i-1)*S+j;
        phi = Find_Eigenfunctions_by_LUT({InputSpace_train{1}(i);InputSpace_train{2}(j)},E);
        tmp = phi' * phi;
        alpha(n,:) = tmp(:);
        beta(n,:) = phi' * y_s(n);
    end
end

% transition matrix
tran = eye(S*S) - Adj;

% Start iteration
f_E = cell(S*S,1);
Pi_E = cell(S*S,1);
m = 1;
while m <= M
    for n = 1:S*S        
        if mod(m,interval) == 0
            if n == SHOW
                f_E{n} = PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(n,:)';
                Pi_E{n} = K - PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * reshape(alpha(n,:),E,E)*LAMBDA*PHI';                
            end       
        end
    end
    % average consensus
    alpha = alpha - gamma*tran*alpha;
    beta = beta - gamma*tran*beta;
    % plot
    if mod(m,interval) == 0
        figure(m+10)
        pcolor(reshape(f_E{SHOW},a,a));        
        colorbar
        figure(m+100)
        pcolor(reshape(diag(Pi_E{SHOW}),a,a));        
        colorbar
    end
    m = m+1;
end