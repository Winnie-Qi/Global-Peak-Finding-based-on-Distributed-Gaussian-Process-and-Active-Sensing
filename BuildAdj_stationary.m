function [distances, tran] = BuildAdj_stationary(locations,r)
% BuildAdj - To compute the adjacency matrix for each agent and plot
% Inputs:
%   locations - (n, 2)matrix, S is the number of sensors/agents, (x,y)
%   coordinates each row
%   r - The range of communication with surrounding sensors
% Output:
%   Adj - The adjacency matrix of the communication network
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
plot(G,'XData',locations(:,2),'YData',locations(:,1));
xlim([-4 4]);
ylim([-4 4]); 
pause(0.01);
hold on
% title('Communication Network');
