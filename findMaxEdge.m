function [goto_x,goto_y] = findMaxEdge(MovingAgent,f,Exploration)

global InputSpace_LUT;
a = size(f,1);
[~, index_x] = min(abs(MovingAgent(1) - InputSpace_LUT{1})); 
[~, index_y] = min(abs(MovingAgent(2) - InputSpace_LUT{2}));

if Exploration
    L = 10; % side length 21    
    while true        
        % Calculates the start and end positions of the window
        start_x = max(index_x-L, 1);
        end_x = min(index_x+L, a);
        start_y = max(index_y-L, 1);
        end_y = min(index_y+L, a);
        % Find the maximum value in the window
        window = f(start_x:end_x, start_y:end_y);
        [max_value, max_index] = max(window(:));
        [max_index_x, max_index_y] = ind2sub(size(window), max_index);
%         goto_x = start_x + max_index_x - 1;
%         goto_y = start_y + max_index_y - 1;
        goto_x = L - max_index_x + 1;
        goto_y = L - max_index_y + 1;
        if max_value > 1
            break
        end
        L = L+10;
    end
else
end
        

end

