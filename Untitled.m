clc;
% close all;
clear all;
num_test_points = 50;
num_train_points = 10;
E = 5; % Dimensions of KL Kernel Expansion
sigma_v = 0.1; % Noise variance of Gaussian process

[x_t, y_t] = meshgrid(linspace(-4, 4, num_test_points));
x_test = [x_t(:), y_t(:)]; % test point coordinates
[x, y] = meshgrid(linspace(-4, -1, num_train_points));
x_train = [x(:), y(:)]; % training point coordinates
% y_train = f(x, y) + sigma_v*randn(num_train_points,num_train_points);
y_train = f(x, y);
y_train = y_train(:); % Generate the height data of the training points

% figure(1)
% pcolor(reshape(y_train,num_train_points,num_train_points));
sigma = 1;
cov_matrix = exp(-pdist2(x_test, x_train).^2 / (2 * sigma^2)); % (2500,64)
cov_matrix_train = exp(-pdist2(x_train, x_train).^2 / (2 * sigma^2)); % (64,64)
y_2 = cov_matrix*pinv(cov_matrix_train+sigma_v^2*num_train_points)*y_train;
cov_matrix_test = exp(-pdist2(x_test, x_test).^2 / (2 * sigma^2));
% figure(2)
% pcolor(cov_matrix)
% shading flat
% figure(3)
% pcolor(pinv(cov_matrix_train+sigma_v^2*num_train_points))
figure(5)
pcolor(x_t,y_t,reshape(y_2,num_test_points,num_test_points));
colorbar
% var_test = cov_matrix_test - cov_matrix * (cov_matrix_train \ cov_matrix');
% figure(1)
% mesh(var_test)
% shading flat
% colorbar
% figure(2)
% mesh(reshape(diag(var_test),num_test_points,num_test_points));
% colorbar

% 
% clc;
% close all;
% clear all;
% 
% % The parameters here can be modified by the users
% InputSpace = {linspace(-4, 4, 100); linspace(-4, 4, 100)};
% PercentageOfVarianceToBeCaptured = 0.96;
% v=0.05;
% l=0.05;
% 
% n_input = numel(InputSpace{1});
% kernels = zeros(n_input);
% 
% % compute the 1D kernel
% for i = 1:n_input
%     for j = 1:n_input
%         i_value = InputSpace{1}(i);
%         j_value = InputSpace{1}(j);
%         kernels(i,j) = exp(-0.5*((i_value-j_value)^2)/v)*l;
%     end
% end
% 
% % compute the eigenfunctions and eigenvalues of the 1D kernel
% [Eigenfunctions_1D, Eigenvalues_1D] = eig(kernels);
% Eigenvalues_1D = abs( Eigenvalues_1D );
% 
% % sort and select the first E eigenfunctions and eigenvalues of the 1D kernel
% [SortedEigenvalues, SortingIndexes] = sort( diag(Eigenvalues_1D), 'descend' );
% SortedEigenfunctions = zeros(size(Eigenfunctions_1D));
% for e = 1:numel(SortingIndexes)    
%     SortedEigenfunctions(:,e) = Eigenfunctions_1D(:, SortingIndexes(e));
% end
% [SortedEigenfunctions_1D,SortedEigenvalues_1D,E] = Select_E(SortedEigenfunctions,SortedEigenvalues,1,PercentageOfVarianceToBeCaptured);
% SortedEigenfunctions_1D = Orthogonalize(SortedEigenfunctions_1D,n_input,1); % Orthonormalize the eigenfunctions
% 
% % combine 1D results to get the 2D results
% [OriginalFirstIndex,OriginalSecondIndex,AllTheCouples]=deal(zeros(E*E, 1));
% 
% for e = 1:E*E
%     iAIndex = ceil( e / E );
%     iBIndex = e - ( (iAIndex - 1) * E );    
%     AllTheCouples(e) = SortedEigenvalues_1D( iAIndex )*SortedEigenvalues_1D( iBIndex );
% end
% 
% [Eigenvalues, aiIndexes] = sort( AllTheCouples, 'descend' );
% 
% for e = 1:E*E
%     OriginalFirstIndex(e)	= ceil( aiIndexes(e)/E);
%     OriginalSecondIndex(e)	= aiIndexes(e) - ( (OriginalFirstIndex(e) - 1) *E);    
% end
% 
% Eigenfunctions = zeros(n_input,n_input,E);
% for e = 1:E*E
%     EigenfunctionIndexInA = OriginalFirstIndex(e);
% 	EigenfunctionIndexInB = OriginalSecondIndex(e);
%     for a = 1:n_input
%         for b = 1:n_input
%             Eigenfunctions(a,b,e) = SortedEigenfunctions_1D( a, EigenfunctionIndexInA )...
%                 *SortedEigenfunctions_1D( b, EigenfunctionIndexInB );
%         end        
%     end
% end
% 
% [Eigenfunctions, Eigenvalues,E] = Select_E(Eigenfunctions,Eigenvalues,2,PercentageOfVarianceToBeCaptured);
% OrthonormalEigenfunctions = Orthogonalize(Eigenfunctions,n_input,2); % Orthonormalize the eigenfunctions
% 
% % Create results directory
% mkdir ./mat
% 
% % 下面是临时的
% S = 10; % The number of stationary sensors is S*S，corresponds to traing points density
% r = 1.8; % The range of communication with surrounding sensors
% a = 100; % test points density
% E = 100;
% sigma = 0.0001;
% M = 100; % iteration times
% SHOW = 6; % show the result of the Sth sensor
% gamma = 0.02;
% InputSpace = {linspace(-4, 4, a); linspace(-4, 4, a)};
% % ShowEnvironment3D(InputSpace);
% ShowTopDownView(InputSpace);
% 
% % Compute the eigenfunctions and eigenvalues of Inputspace
% PHI = OrthonormalEigenfunctions(:,:,1:E);
% PHI = reshape(PHI,[],E);
% LAMBDA = diag(Eigenvalues(1:E));
% 
% [temp_x, temp_y] = meshgrid(linspace(-4, 4, S));
% StationarySensors = zeros(S*S, 2);
% StationarySensors(:, 1) = temp_x(:);
% StationarySensors(:, 2) = temp_y(:);
% Adj = BuildAdj(StationarySensors, r); % Adjacency matrix
% y_s = f(temp_x, temp_y) + sigma*randn(S,S);
% y_s = y_s(:);
% G = OrthonormalEigenfunctions(1:10:100,1:10:100,1:E);
% G = reshape(G,[],E);
% 
% f_E = PHI * pinv(G'*G/S^2+sigma^2/S^2*pinv(LAMBDA)) * G'/S^2*y_s; 
% figure(3)
% pcolor(reshape(f_E,a,a));
% colorbar