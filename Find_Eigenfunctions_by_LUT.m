function EigFuns = Find_Eigenfunctions_by_LUT(InputSpace_1,E)

global Eig_LUT InputSpace_LUT;
n_input_1 = numel(InputSpace_1{1});
n_input_2 = numel(InputSpace_1{2});
EigFuns = zeros(n_input_1,n_input_2,E);
for i = 1:n_input_1
    % Get the nearest neighbor index of the current coordinate in InputSpace
    [~, index_i] = min(abs(InputSpace_1{1}(i) - InputSpace_LUT{1}));
    for j = 1:n_input_2        
        [~, index_j] = min(abs(InputSpace_1{2}(j) - InputSpace_LUT{2}));        
        % Copy the corresponding Eigenfunctions value to PHI
        EigFuns(i, j, :) = Eig_LUT(index_i,index_j,1:E);
    end
end
EigFuns = reshape(EigFuns,[],E);
end

% global Eig_LUT InputSpace_LUT;
% n_input_1 = numel(InputSpace_1{1});
% n_input_2 = numel(InputSpace_1{2});
% EigFuns = zeros(n_input_1,n_input_2,E);
% for i = 1:n_input_1
%     % Get the nearest neighbor index of the current coordinate in InputSpace
%     [~, index_i] = min(abs(InputSpace_1{1}(i) - InputSpace_LUT{1}));
%     for j = 1:n_input_2        
%         [~, index_j] = min(abs(InputSpace_1{2}(j) - InputSpace_LUT{2}));        
%         % Copy the corresponding Eigenfunctions value to PHI
%         EigFuns(i, j, :) = Eig_LUT(index_i,index_j,1:E);
%     end
% end
% EigFuns = reshape(EigFuns,[],E);
% end
