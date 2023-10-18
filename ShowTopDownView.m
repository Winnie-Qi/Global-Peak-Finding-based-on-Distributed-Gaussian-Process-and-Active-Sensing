function ShowTopDownView(InputSpace)
% True environment from the top down view

[x,y] = meshgrid(InputSpace{1},InputSpace{2});
z =  f(x',y');
h2 = pcolor(x,y,z);
h2.EdgeColor = 'none';
title('Top down view');
colorbar
pause(0.01);
hold on

