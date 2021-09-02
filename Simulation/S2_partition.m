% scenario 2 
% Partition algorithm

function [partition_result, number_of_subgroups, middle_point] = S2_partition(X, Y, UAVradius)
    n_of_nodes = length(X);
    number_of_subgroups = n_of_nodes;
    
    % Distance matrix
    pairwiseD = zeros(n_of_nodes, n_of_nodes);
    for i = 1: n_of_nodes
        for j = 1: n_of_nodes
            pairwiseD(i, j) = ( ( X(i) - X(j) )^2 +  ( Y(i) - Y(j) )^2 )^(1/2); 
        end
    end
    
    % calculate distance between each pair,    use a m*3 matrix to store the distance between each pair [index1, index2, distance]
    Distance = zeros( (n_of_nodes^2 - n_of_nodes)/2  , 3);
    tempnum = 1;
    for i = 1: n_of_nodes-1
        for j = i+1: n_of_nodes
            distance = ( ( X(i) - X(j) )^2 +  ( Y(i) - Y(j) )^2 )^(1/2);
            Distance( tempnum, :) = [i, j,  distance]; 
            tempnum = tempnum +1 ;
        end
    end

    Distance = sortrows(Distance, 3);   % sort in increasing order
    
    % grow spanning tree,     use a 2*n matrix to perform the set operation (the set state),      [  ID ;  set ID ]
    set_node = [  1:length(X) ;  1:length(X)  ]; 
    set_element = zeros(length(X), length(X) + 1);     % store the element in set i, the last column is the number of element 
    
    % initialization
    set_element(:,1) = 1:length(X);                    % element i is in set i
    set_element(:, length(X) + 1) = ones(length(X) , 1);
    
    for i = 1:length(Distance(:, 1))
        
        edge = Distance(i, :);    % [index1, index2, distance]
        
        if edge(3) > UAVradius*(3^(1/2))   % cannot fulfill, and since the edge is increasing, we can break here
            break;
        end
        
        % check if this forms a cycle        
        setIndex1 = set_node(  2, edge(1)   );  % node index = edge(1) and edge(2)
        setIndex2 = set_node(  2,  edge(2)   );
        if setIndex1 ~= setIndex2               % not in the same subtree 
                                               
            merge_or_not = 1;              % check if new diameter is within the range
            for ti = 1: set_element(setIndex1, length(X)+1)
                tempN1 = set_element(setIndex1,ti);
                if tempN1 == 0   % no more element in set 1
                    break;
                end
                for tj = 1:1: set_element(setIndex2, length(X)+1)
                    tempN2 = set_element(setIndex2,tj);
                    if tempN2 == 0   % no more element in set 2, break;
                        break;
                    end 
                    
                    if pairwiseD(tempN2, tempN1) > UAVradius*(3^(1/2)) % diamete exceeds, cannot merge
                        merge_or_not = 0;
                        break;
                    end
                end
                if merge_or_not == 0
                    break;
                end
            end            % check if new diameter is within the range
                        
            if merge_or_not == 0       % out of diameter, continue trying next one
                continue;
            else                      % number of groups --    
                number_of_subgroups = number_of_subgroups -1;
            end
            
            % merge two sets (change setindex2)  
            % 
            for k = 1:set_element(setIndex2, length(X)+1 ) 
                temp_node_index = set_element(setIndex2, k);
                set_element(setIndex2, k) = 0;
                set_node(2, temp_node_index ) =  setIndex1;
                set_element(setIndex1, set_element(setIndex1, length(X)+1 ) +  k) = temp_node_index ;
            end
            set_element(setIndex1, length(X)+1 )  = set_element(setIndex2, length(X)+1 ) + set_element(setIndex1, length(X)+1 );
            set_element(setIndex2, length(X)+1 )  = 0;
            
        else  % this edge forms a cycle, thus, skip this edge
            continue;    
        end
    end
    
    % reorder the set index
    % set_node 
    unusedindex =1;
    %set_node(2,:) = set_node(2,:) -100000;   % make sure the label won't be the same as 1,2,3......
    
    hash_table = -1*ones(1, 10000);   % -1 means no hit history, o.w., gives an unusedindex
    
    for iind = 1:length(X)
        cur_subg_index = set_node(2, iind);
        if hash_table( cur_subg_index ) == -1
            hash_table(cur_subg_index) = unusedindex;
            set_node(2, iind) = unusedindex;
            unusedindex = unusedindex + 1;
        else
            set_node(2, iind) = hash_table(cur_subg_index);
        end
    end

    
    % --------------------------------------------------- calculate the middle point of each subgroup for tree partition -----------------------------------
    % index stard from 1
    num_subgroups_temp = number_of_subgroups;
    partition_resilt_temp = set_node;
    middle_point = -1000*ones(2, num_subgroups_temp);           
                                                             
    for subgroup_index = 1:num_subgroups_temp 
        index_of_devices_sbg = find(partition_resilt_temp(2,:) == subgroup_index);    % return the indexes of devices in this group
        
        if isempty(index_of_devices_sbg)
            fprintf('ERROR, should not contains empty set with index < num_subgroups_temp\n');
            pause
            %break;    % since the subgroup index is sorted and reordered, stop when reaching the first empty subgroups
        end
        mid_x = mean( X(index_of_devices_sbg) );
        mid_y = mean( Y(index_of_devices_sbg) );
        middle_point(:, subgroup_index) = [mid_x; mid_y];
    end
    
    
    
    partition_result = set_node;
end


%{

    for ind = 1:length(X)
        temp_sbg_index = set_node(2, ind);
        if temp_sbg_index > 0 % already reordered
            continue;
        end
        % find all the devices in the same group
        index_of_devices_sbg = find( set_node(2,:) == temp_sbg_index );
        
        if isempty(index_of_devices_sbg) % this index is not used
            fprintf('ERROR!!!');
            pause;
        end
            
        % change index for these devices
        for ii = 1:length(index_of_devices_sbg)
            tind = index_of_devices_sbg(ii);
            set_node(2, tind) = unusedindex;
        end
        unusedindex = unusedindex +1;
        if unusedindex > number_of_subgroups
            break;
        end
    end

%}




%{

 Large_set_index = max(setIndex1, setIndex2);
            small_set_index = min(setIndex1, setIndex2);
            for k = 1:set_element(Large_set_index, length(X)+1 ) 
                temp_node_index = set_element(Large_set_index, k);
                set_element(Large_set_index, k) = 0;
                set_node(2, temp_node_index ) =  small_set_index;
                
                set_element(small_set_index, set_element(small_set_index, length(X)+1 ) +  k) = temp_node_index ;
            end
            set_element(small_set_index, length(X)+1 )  = set_element(Large_set_index, length(X)+1 ) + set_element(small_set_index, length(X)+1 );
            set_element(Large_set_index, length(X)+1 )  = 0;
%}