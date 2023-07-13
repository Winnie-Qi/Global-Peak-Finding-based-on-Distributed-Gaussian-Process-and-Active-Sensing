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
E = 6;
sigma = 0.0001;
M = 100; % iteration times
SHOW = 6; % show the result of the Sth sensor
gamma = 0.02;

%%

% Show the environment
InputSpace = {linspace(-4, 4, a); linspace(-4, 4, a)};
% ShowEnvironment3D(InputSpace);
ShowTopDownView(InputSpace);

% Compute the eigenfunctions and eigenvalues of Inputspace
PHI = zeros(a^2,E);
t_vals = 1:a;
s_vals = 1:a;
PHI(:,1) = reshape(2*sin(1/2*pi*InputSpace{1}(s_vals))' * sin(1/2*pi*InputSpace{2}(t_vals)),1,[]);
PHI(:,2) = reshape(2*sin(1/2*pi*InputSpace{1}(s_vals))' * sin(3/2*pi*InputSpace{2}(t_vals)),1,[]);
PHI(:,3) = reshape(2*sin(3/2*pi*InputSpace{1}(s_vals))' * sin(1/2*pi*InputSpace{2}(t_vals)),1,[]);
PHI(:,4) = reshape(2*sin(5/2*pi*InputSpace{1}(s_vals))' * sin(1/2*pi*InputSpace{2}(t_vals)),1,[]);
PHI(:,5) = reshape(2*sin(1/2*pi*InputSpace{1}(s_vals))' * sin(5/2*pi*InputSpace{2}(t_vals)),1,[]);
PHI(:,6) = reshape(2*sin(1/2*pi*InputSpace{1}(s_vals))' * sin(7/2*pi*InputSpace{2}(t_vals)),1,[]);
% for i = 1:1
%     for j = 1:E
%         phi = 2 * sin((i-0.5)*pi*InputSpace{1})' * sin((j-0.5)*pi*InputSpace{2});
%         PHI(:, (i-1)*E+j) = reshape(phi,1,[]);
%     end
% end

LAMBDA = zeros(E,E);
LAMBDA(1,1) = 16/(1/2^2*1/2^2*pi^4);
LAMBDA(2,2) = 16/(1/2^2*3/2^2*pi^4);
LAMBDA(3,3) = 16/(3/2^2*1/2^2*pi^4);
LAMBDA(4,4) = 16/(5/2^2*1/2^2*pi^4);
LAMBDA(5,5) = 16/(1/2^2*5/2^2*pi^4);
LAMBDA(6,6) = 16/(1/2^2*7/2^2*pi^4);
% LAMBDA = diag(1 ./ (((2*(1:E)-1).^2*pi^2).^2));

% Initialize the stationary sensors
[temp_x, temp_y] = meshgrid(linspace(-3.5, 3.5, S));
StationarySensors = zeros(S*S, 2);
StationarySensors(:, 1) = temp_x(:);
StationarySensors(:, 2) = temp_y(:);
Adj = BuildAdj(StationarySensors, r); % Adjacency matrix
y_s = f(temp_x, temp_y) + sigma*randn(S,S);
y_s = y_s(:);

G = zeros(S*S,E);
% for i = 1:1
%     for j = 1:E
%         phi = 2 * sin((i-0.5)*pi*StationarySensors(:,1)) .* sin((j-0.5)*pi*StationarySensors(:,2));
%         lambda = 16 / (((2*i-1)^2*pi^2) * ((2*j-1)^2*pi^2));
%         G(:, (i-1)*E+j) = phi * sqrt(lambda);
%     end
% end
for n = 1:S*S
    G(n,1) = 2*sin(1/2*pi*StationarySensors(n,1)) * sin(1/2*pi*StationarySensors(n,2));
    G(n,2) = 2*sin(1/2*pi*StationarySensors(n,1)) * sin(3/2*pi*StationarySensors(n,2));
    G(n,3) = 2*sin(3/2*pi*StationarySensors(n,1)) * sin(1/2*pi*StationarySensors(n,2));
    G(n,4) = 2*sin(5/2*pi*StationarySensors(n,1)) * sin(1/2*pi*StationarySensors(n,2));
    G(n,5) = 2*sin(1/2*pi*StationarySensors(n,1)) * sin(5/2*pi*StationarySensors(n,2));
    G(n,6) = 2*sin(1/2*pi*StationarySensors(n,1)) * sin(7/2*pi*StationarySensors(n,2));   
end
f_E = PHI * pinv(G'*G/S^2+sigma^2/S^2*pinv(LAMBDA)) * G'/S^2*y_s; 
figure(3)
pcolor(reshape(f_E,a,a));

% Intialize alpha and beta
alpha = zeros(S*S,E*E);
beta = zeros(S*S,E);
PHI = zeros(E,1);
for n = 1:S*S    
    for e = 1:E
        PHI(e) = 2*sin((e-1/2)*pi*StationarySensors(n,1)) * 2*sin((e-1/2)*pi*StationarySensors(n,2));
    end
    tmp = PHI * PHI';
    alpha(n,:) = tmp(:);
    beta(n,:) = PHI * y_s(n);
end

% transition matrix
tran = eye(S*S) - Adj;

% Start iteration
f_E = cell(S*S,1);
m = 1;
while m <= M
    for n = 1:S*S
        f_E{n} = PHI * pinv(reshape(alpha(n,:),E,E)+sigma^2/S^2*pinv(LAMBDA)) * beta(n,:)'; 
%         f_E{n} = PHI * eye(E,E)/(reshape(alpha(n,:),E,E)+sigma^2/S^2*eye(E,E)/LAMBDA) * beta(n,:)'; 
        alpha = alpha + gamma*(tran*alpha);
        beta = beta + gamma*(tran*beta);
    end
    if mod(m,20) == 0
        figure(m)
        pcolor(reshape(f_E{6},a,a));
    end
    m = m+1;
end