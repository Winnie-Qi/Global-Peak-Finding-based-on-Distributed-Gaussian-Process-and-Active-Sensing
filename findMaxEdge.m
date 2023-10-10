function [goto_x,goto_y,ExplorationFlag,reachGoal] = findMaxEdge(MovingAgent,f,Exploration)

global InputSpace_LUT;
ExplorationFlag = 1;
reachGoal = 0;
m = 1;
M = 9;
a = size(f,1);
[~, index_x] = min(abs(MovingAgent(1) - InputSpace_LUT{1})); 
[~, index_y] = min(abs(MovingAgent(2) - InputSpace_LUT{2}));

if Exploration
    L = 2; % side length 21    
    while true        
        % Calculates the start and end positions of the window
        start_x = max(index_x-L, 1);
        end_x = min(index_x+L, a);
        start_y = max(index_y-L, 1);
        end_y = min(index_y+L, a);
        % Find the maximum value in the window
        window = f(start_x:end_x, start_y:end_y);
        [max_value, max_index] = max(window(:));
%         disp(max_value)
        [max_index_x, max_index_y] = ind2sub(size(window), max_index);
%         goto_x = start_x + max_index_x - 1;
%         goto_y = start_y + max_index_y - 1;
        goto_x = L - max_index_x + 1;
        goto_y = L - max_index_y + 1;
        if max_value > M
            s = sqrt(goto_x^2 + goto_y^2);
            if s == 0
                goto_x = 1; goto_y = 1;
                m = 1;
            else
                goto_x = goto_x/s * m;
                goto_y = goto_y/s * m;
            end
            break
        else
            if L > 7
                goto_x=0; goto_y=0;
                ExplorationFlag = 0;
                break
            end
        end
        L = L+1;
        m = m + 1;
        M = M - 0.5;        
    end
else
    [~, max_index] = max(f(:));
    [max_index_x, max_index_y] = ind2sub(size(f), max_index);
%     disp(max_index_x)
%     disp(max_index_y)
%     disp(max_value)
    goto_x = max_index_x - index_x + 1;
    goto_y = max_index_y - index_y + 1;
    s = sqrt(goto_x^2 + goto_y^2);
    if s == 0
        goto_x = 1; goto_y = 1;        
    else
        goto_x = goto_x/s;
        goto_y = goto_y/s;
    end
    if s < 12
        reachGoal = 1;
    end
        
end        

end