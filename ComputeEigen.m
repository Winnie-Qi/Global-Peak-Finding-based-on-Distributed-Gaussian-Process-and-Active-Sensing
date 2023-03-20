function [OrthonormalEigenfunctions, Eigenvalues] = ComputeEigen(InputSpace,v,l)
% To simpify the computation,we assume that InputSpace is symmetric based on the origin
if nargin < 2
    v=0.05;
    l=0.05;
end
n_input = numel(InputSpace{1});
kernels = zeros(n_input);

% compute the 1D kernel
for i = 1:n_input
    for j = 1:n_input
        i_value = InputSpace{1}(i);
        j_value = InputSpace{1}(j);
        kernels(i,j) = exp(-0.5*((i_value-j_value)^2)/v)*l;
    end
end

% compute the eigenfunctions and eigenvalues of the 1D kernel
[Eigenfunctions_1D, Eigenvalues_1D] = eig(kernels);
Eigenvalues_1D = abs( Eigenvalues_1D );

% sort and select the first E eigenfunctions and eigenvalues of the 1D kernel
[SortedEigenvalues, SortingIndexes] = sort( diag(Eigenvalues_1D), 'descend' );
SortedEigenfunctions = zeros(size(Eigenfunctions_1D));
for e = 1:numel(SortingIndexes)    
    SortedEigenfunctions(:,e) = Eigenfunctions_1D(:, SortingIndexes(e));
end
[SortedEigenfunctions_1D,SortedEigenvalues_1D,E] = Select_E(SortedEigenfunctions,SortedEigenvalues,1);
% [SortedEigenvalues_1D, SortingIndexes] = sort( diag(Eigenvalues_1D), 'descend' );
% PropotionSorting = cumsum(SortedEigenvalues_1D)./sum(SortedEigenvalues_1D);
% E = find((PropotionSorting < 0.995),1,'last');
% fprintf('E is: %d in one dimension.\n', E*E);
% SortedEigenfunctions_1D = zeros(numel(InputSpace{1}),E);
% for e = 1:E
%     SortedEigenfunctions_1D(:,e) = Eigenfunctions_1D(:, SortingIndexes(e));
% end
SortedEigenfunctions_1D = Orthogonalize(SortedEigenfunctions_1D,n_input,1); % Orthonormalize the eigenfunctions
% SortedEigenvalues_1D = SortedEigenvalues_1D(1:E);

% combine 1D results to get the 2D results
[OriginalFirstIndex,OriginalSecondIndex,AllTheCouples]=deal(zeros(E*E, 1));

for e = 1:E*E
    iAIndex = ceil( e / E );
    iBIndex = e - ( (iAIndex - 1) * E );    
    AllTheCouples(e) = SortedEigenvalues_1D( iAIndex )*SortedEigenvalues_1D( iBIndex );
end

[Eigenvalues, aiIndexes] = sort( AllTheCouples, 'descend' );

for e = 1:E*E
    OriginalFirstIndex(e)	= ceil( aiIndexes(e)/E);
    OriginalSecondIndex(e)	= aiIndexes(e) - ( (OriginalFirstIndex(e) - 1) *E);    
end

Eigenfunctions = zeros(n_input,n_input,E);
for e = 1:E*E
    EigenfunctionIndexInA = OriginalFirstIndex(e);
	EigenfunctionIndexInB = OriginalSecondIndex(e);
    for a = 1:n_input
        for b = 1:n_input
            Eigenfunctions(a,b,e) = SortedEigenfunctions_1D( a, EigenfunctionIndexInA )...
                *SortedEigenfunctions_1D( b, EigenfunctionIndexInB );
        end        
    end
end

[Eigenfunctions, Eigenvalues,E] = Select_E(Eigenfunctions,Eigenvalues,2);
OrthonormalEigenfunctions = Orthogonalize(Eigenfunctions,n_input,2); % Orthonormalize the eigenfunctions