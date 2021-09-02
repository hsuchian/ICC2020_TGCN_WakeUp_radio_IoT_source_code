% even partition

function [partition_result, number_of_subgroups, middle_point] = S2_EvenPartition_square(X, Y, UAVradius, max_range)
    
    n_of_nodes = length(X);
    % number_of_subgroups = n_of_nodes;
    
    set_node = [  1:length(X) ;  -1*ones(1, n_of_nodes) ]; 
    hash_table_1 = zeros(1, 10000);   % -1 means no hit history, o.w., gives an unusedindex
    
    width = 2^(1/2)*UAVradius;
    for i = 1:n_of_nodes
        
        x_ind = ceil( X(i)/width); 
        y_ind = ceil( Y(i)/width); 
        set_node(2, i) = x_ind + (y_ind-1)*ceil(max_range/width); 
        hash_table_1(set_node(2, i)) = 1;
    end
    
    num_subgroups = sum(hash_table_1);
    
    % calculate the middle points  
    middle_point = -1*ones(2, num_subgroups);   %[x;y]
    
    % reorder the set index
    % set_node 
    unusedindex =1;
    %set_node(2,:) = set_node(2,:) -100000;   % make sure the label won't be the same as 1,2,3......
    
    set_node_old = set_node;
    
    hash_table = -1*ones(1, 10000);   % -1 means no hit history, o.w., gives an unusedindex
    sbg_bun_side = ceil(max_range/width);
    for iind = 1:length(X)
        cur_subg_index = set_node(2, iind);
        if hash_table( cur_subg_index ) == -1
            % first hit, calculate the middle point of this group according to the real group index
            % the mid point is stored according to the new assigned index
            mid_x = mod(cur_subg_index, sbg_bun_side );
            if mid_x ==0
                 mid_x = sbg_bun_side*width - width/2; 
            else
                mid_x = mid_x*width - width/2; 
            end
            mid_y = ceil( cur_subg_index/ sbg_bun_side );
            mid_y = mid_y*width - width/2; 
            
            % update the index
            hash_table(cur_subg_index) = unusedindex;
            set_node(2, iind) = unusedindex;
            unusedindex = unusedindex + 1;
            
            middle_point(1, set_node(2, iind)) = mid_x;
            middle_point(2, set_node(2, iind)) = mid_y;
        else
            set_node(2, iind) = hash_table(cur_subg_index);
        end
    end

    number_of_subgroups = unusedindex-1;
    partition_result = set_node;
end 