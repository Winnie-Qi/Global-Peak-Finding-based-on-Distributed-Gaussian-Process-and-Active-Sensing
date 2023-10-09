function [distances, tran] = BuildAdj(locations,r,InputSpace_test,S)
% BuildAdj - To compute the adjacency matrix for each agent and plot
% Inputs:
%   locations - (n, 2)matrix, S is the number of sensors/agents, (x,y)
%   coordinates each row
%   r - The range of communication with surrounding sensors
% Output:
%   Adj - The adjacency matrix of the communication network
subplot(4,S,[4,5,6,4+S,5+S,6+S]);
ShowTopDownView(InputSpace_test);
n = size(locations,1); % number of agents/sensors
Adj = zeros(n);% adjacency matrix with all initial values of 0
distances = zeros(n);
pre_tran = zeros(1,n);
for i =1:n-1
    for j = i+1:n
        d = locations(j,:) - locations(i,:); % The distance between every two points
        d = sqrt(sum(d.*d));        
        if d <= r
            Adj(i,j)=1;
            Adj(j,i)=1;
            pre_tran(1,i)=1;
            if d <= 0.6
                distances(i,j)=1;
                distances(j,i)=1;
            end
        end
    end
end
tran = diag(pre_tran) - Adj;
G = digraph(Adj);
% figure(3)
plot(G,'XData',locations(:,1),'YData',locations(:,2));
xlim([-4 4]);
ylim([-4 4]); 
% title('Communication Network');
z = ones(S,1) * 15;
subplot(4,S,[1,2,3,1+S,2+S,3+S]);
ShowEnvironment3D(InputSpace_test);
plot(G,'XData',locations(:,1),'YData',locations(:,2),'ZData',z);
hold off