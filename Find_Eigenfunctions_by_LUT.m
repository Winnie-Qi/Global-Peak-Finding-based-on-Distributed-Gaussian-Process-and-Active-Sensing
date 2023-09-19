function EigFuns = Find_Eigenfunctions_by_LUT(InputSpace_1,E)

global Eig_LUT InputSpace_LUT;

n_input = size(InputSpace_1,1);
EigFuns = zeros(n_input,E);
for i = 1:n_input
    [~, index_1] = min(abs(InputSpace_1(i,1) - InputSpace_LUT{1})); 
    [~, index_2] = min(abs(InputSpace_1(i,2) - InputSpace_LUT{2})); 
    EigFuns(i, :) = Eig_LUT(index_1,index_2,1:E);   
end