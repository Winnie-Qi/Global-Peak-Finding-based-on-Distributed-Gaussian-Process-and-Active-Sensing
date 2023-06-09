function Adj = BuildAdj(locations,r)
% BuildAdj - To compute the adjacency matrix for each agent and plot
% Inputs:
%   locations - (n, 2)matrix, S is the number of sensors/agents, (x,y)
%   coordinates each row
%   r - The range of communication with surrounding sensors
% Output:
%   Adj - The adjacency matrix of the communication network

n = size(locations,1); % number of agents/sensors
Adj = zeros(n); % adjacency matrix with all initial values of 0
for i =1:n-1
    for j = i+1:n
        d = locations(j,:) - locations(i,:); % The distance between every two points
        d = sqrt(sum(d.*d));
        if d <= r
            Adj(i,j)=1;
            Adj(j,i)=1;
        end
    end
end
G = digraph(Adj);
figure
plot(G,'XData',locations(:,1),'YData',locations(:,2));
xlim([-4 4]);
ylim([-4 4]);
title('Communication Network');