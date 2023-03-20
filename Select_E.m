function [SortedEigenfunctions,SortedEigenvalues,E] = Select_E(Eigenfunctions,Eigenvalues,d)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
PropotionSorting = cumsum(Eigenvalues)./sum(Eigenvalues);
E = find((PropotionSorting < 0.9),1,'last');
fprintf('E is %d in %d dimension.\n', E,d);
switch d
    case 1        
        SortedEigenfunctions = Eigenfunctions(:,1:E);
    case 2
        SortedEigenfunctions = Eigenfunctions(:,:,1:E);
end
SortedEigenvalues = Eigenvalues(1:E);
end

