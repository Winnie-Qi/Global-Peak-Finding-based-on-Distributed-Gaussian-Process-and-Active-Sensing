function OrthogonalEigenfunctions = Orthogonalize(SortedEigenfunctions, n_input, d)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
OrthogonalEigenfunctions = zeros( size(SortedEigenfunctions) );
switch d
    case 1
        E = size(SortedEigenfunctions,2);
        for e = 1:E
            CurrentEigenfunctionL2MuNorm = sum(SortedEigenfunctions(:,e).^2)/n_input;
            OrthogonalEigenfunctions(:, e) = SortedEigenfunctions(:,e)./sqrt(CurrentEigenfunctionL2MuNorm);
        end
    case 2
        E = size(SortedEigenfunctions,3);
        for e = 1:E
            CurrentEigenfunctionL2MuNorm = sum(sum(SortedEigenfunctions(:,:,e).^2./n_input^2));
			OrthogonalEigenfunctions(:,:,e) = SortedEigenfunctions(:,:,e)./sqrt(CurrentEigenfunctionL2MuNorm);
        end
end