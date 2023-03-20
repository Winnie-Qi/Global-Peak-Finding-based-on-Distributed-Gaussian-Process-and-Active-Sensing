function ShowTopDownView(InputSpace)
% True environment from the top down view
figure
[x,y] = meshgrid(InputSpace{1},InputSpace{2});
z =  f(x,y);
p = pcolor(x,y,z);
p.EdgeColor = 'none';
title('3D environmet from top down view');
colorbar

end

