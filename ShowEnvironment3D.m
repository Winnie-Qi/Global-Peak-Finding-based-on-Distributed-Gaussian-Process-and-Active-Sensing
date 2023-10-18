function ShowEnvironment3D(InputSpace)

[x,y] = meshgrid(InputSpace{1},InputSpace{2});
z =  f(x',y');
mesh(x,y,z);
title('3D Environmet');
xlim([-4 4]);
ylim([-4 4]);
zlim([-12.5 20]);
pause(0.01);
hold on