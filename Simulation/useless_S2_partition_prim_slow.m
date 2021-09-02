% scenario 2 
% Partition algorithm

function [partition_result, number_of_subgroups, middle_point] = S2_partition_prim(X, Y, UAVradius)
    
    int_max = 2147483647;
    n_of_nodes = length(X);
    
    
    % use a 2*n matrix to perform the set operation (the set state),      [  ID ;  set ID ]
    set_node = [  1:length(X) ;  1:length(X)  ]; 
    set_element = zeros(length(X), length(X) + 1);     % store the element in set i, the last column is the number of element 
    
    
    % Distance matrix
    pairwiseD = zeros(n_of_nodes, n_of_nodes);
    for i = 1: n_of_nodes
        for j = 1: n_of_nodes
            pairwiseD(i, j) = ( ( X(i) - X(j) )^2 +  ( Y(i) - Y(j) )^2 )^(1/2); 
        end
    end
    
    
    
    % calculate distance (edge) between each pair,    use a m*3 matrix to store the distance between each pair [index1, index2, distance]
    remaining_edge = zeros( (n_of_nodes^2 - n_of_nodes)/2  , 3);    % the set of edges E
    edge_num = 1;
    for i = 1: n_of_nodes-1
        for j = i+1: n_of_nodes
            edge_length = ( ( X(i) - X(j) )^2 +  ( Y(i) - Y(j) )^2 )^(1/2);
            remaining_edge( edge_num, :) = [i, j,  edge_length]; 
            edge_num = edge_num +1 ;
        end
    end
    % remaining_edge = sortrows(remaining_edge, 3);       % the remaining edge set
    % Distance = sortrows(Distance, 3);   % sort in increasing order
  
    
    % -------------------------------------- new partition algorithm ------------------------------
                     % [ vertex index ; distance to the current subgroup ; in current subgroup or not  ]
    remaining_vertex = [ 1:n_of_nodes ; int_max*ones(1, n_of_nodes)     ; zeros(1, n_of_nodes)  ];   
    
    %[num_edges, col_size_RE] = size(remaining_edge);
    
    %if col_size_RE ~= 3
    %    fprintf('Error size of Distance matrix n (line 46)\n');
    %    pause 
    %elseif num_edges ~= (n_of_nodes^2 - n_of_nodes)/2
    %    fprintf('Error size of Distance matrix m (line 49)\n');
    %    pause 
    %end 
    
    set_count = 0;
    while isempty(remaining_vertex)  ~= 1
        set_count = set_count+1;
        
        %initialize remaining_vertex
        remaining_vertex(2,:) = int_max*ones(1, length(remaining_vertex(1,:)) );
        remaining_vertex(3,:) = zeros(1, length(remaining_vertex(1,:)) );
        
        %fprintf('Set %d\n', set_count); % debug
        if (isempty(remaining_edge)  ~= 1)
            [min_edge_val, min_edge_ind] = min( remaining_edge(:, 3) );
            min_vertex_ind = remaining_edge(min_edge_ind, 1);   % choose one endpoint as starting point
            %fprintf('    min_edge_val = %d, min_edge_ind = %d\n', min_edge_val, min_edge_ind);
        else 
            min_vertex_ind = remaining_vertex(1,1); % only one vertex remain, the edge set (Distance) is empty
            if length(remaining_vertex(1,:)) ~= 1
                fprintf('Error size of remaining_vertex while facing empty edge set\n');
            end
        end
        
        % add min_vertex_ind into subgroup
        [min_vertex_row_index_in_RV, min_vertex_col_ind_in_RV] = find( remaining_vertex(1,:) == min_vertex_ind  );
        
        if isempty(min_vertex_col_ind_in_RV) || min_vertex_row_index_in_RV ~= 1
            fprintf('ERROR when getting index of min_x line 75\n');
            fprintf('%d   %d\n', min_vertex_ind, min_edge_ind);
            pause
        end
        
        remaining_vertex(2, min_vertex_col_ind_in_RV) = 0;   % useless
        remaining_vertex(3, min_vertex_col_ind_in_RV) = 1;   % useless
        
        if (min_vertex_ind ~= remaining_vertex(1, min_vertex_col_ind_in_RV))
            fprintf('min_vertex_ind %d ~= remaining_vertex(1, min_vertex_col_ind_in_RV) %d\n', min_vertex_ind, remaining_vertex(1, min_vertex_col_ind_in_RV));
        end
        
        set_node(2, min_vertex_ind) = set_count;
        set_element(set_count, n_of_nodes+1) = set_element(set_count, n_of_nodes+1) + 1;
        n_nodes_in_set_count = set_element(set_count, n_of_nodes+1);
        set_element(set_count, n_nodes_in_set_count) =  min_vertex_ind;
        
        % remove the chosen vertex from remaining_vertex
        remaining_vertex(:, min_vertex_col_ind_in_RV) = [];
        
        % update the distance table
        for ind = 1:length(remaining_vertex(1,:))
            if (remaining_vertex(3, ind) == 0)                      % only those remaining vertexes that are not added into current subgroup
                if ( remaining_vertex(2, ind) > pairwiseD( remaining_vertex(1, ind), min_vertex_ind) )
                    remaining_vertex(2, ind) = pairwiseD( remaining_vertex(1, ind), min_vertex_ind);
                end
            else
                % error, should not have any vertex in currnt subgroup
                 fprintf('ERROR, %d is in current subgroup %d (line 85)\n', ind, set_count);
            end
        end
        
        % remove the adjacent edges of this vertex
        [tmp_num_deges, dummy] = size(remaining_edge);
        ind = 1;
        while ind <= tmp_num_deges
            if remaining_edge(ind, 1) == min_vertex_ind || remaining_edge(ind, 2) == min_vertex_ind
                remaining_edge(ind,:) = [];
                tmp_num_deges = tmp_num_deges - 1;
            else
                ind =ind + 1;
            end
        end
        
        % starting point of inner while loop
        tmp_remaining_vertex = remaining_vertex;   % copy current vertex for current search
                    
        
        while isempty(tmp_remaining_vertex) ~= 1
            [min_value, current_min_vertex_col_ind_in_tmp_RV] = min(tmp_remaining_vertex(2,:));
            
            current_min_vertex_index = tmp_remaining_vertex(1, current_min_vertex_col_ind_in_tmp_RV);
            
            if min_value ~= tmp_remaining_vertex(2, current_min_vertex_col_ind_in_tmp_RV)
                fprintf('Error in choosing min vertex\n');
                pause
            end
            
            % check diameter of S_ik after adding this vertex
            % check the maximum edge between this vertex and S_ik
            feasible = 1;
            for ind = 1:set_element(set_count,  n_of_nodes+1)
                vertex_ind_in_S_ik = set_element(set_count, ind);
                if  ( pairwiseD( vertex_ind_in_S_ik, current_min_vertex_index) >  UAVradius*(3^(1/2))  ) 
                    feasible = 0;
                    break;
                end
            end
            
            if feasible == 1  % the chosen vertex is feasible,
                % add this vertex in curent subgroup
                tmp_remaining_vertex(:, current_min_vertex_col_ind_in_tmp_RV) = [];
                
                [current_min_vertex_row_ind_in_RV, current_min_vertex_col_ind_in_RV] = find( remaining_vertex(1,:) == current_min_vertex_index  );
                remaining_vertex(:, current_min_vertex_col_ind_in_RV ) = [];
        
                set_node(2, current_min_vertex_index) = set_count;
                set_element(set_count, n_of_nodes+1) = set_element(set_count, n_of_nodes+1) + 1;
                n_nodes_in_set_count = set_element(set_count, n_of_nodes+1);
                set_element(set_count, n_nodes_in_set_count) =  current_min_vertex_index;
                
                % update distance table in tmp_remaining_vertex only
                 for ind = 1:length(tmp_remaining_vertex(1,:))
                    if (tmp_remaining_vertex(3, ind) == 0)                      % only those remaining vertexes that are not added into current subgroup
                        if ( tmp_remaining_vertex(2, ind) > pairwiseD( tmp_remaining_vertex(1, ind),  current_min_vertex_index) )
                             tmp_remaining_vertex(2, ind) = pairwiseD( tmp_remaining_vertex(1, ind),  current_min_vertex_index);
                        end
                    else
                        % error, should not have any vertex in currnt subgroup
                        fprintf('ERROR, %d is in current subgroup %d (line 165)\n', ind, set_count);
                    end
                 end  
                 
                 % remove all the edges that is adjacent to current vertex
                 [tmp_num_deges, dummy] = size(remaining_edge);
                 ind = 1;
                 while ind <= tmp_num_deges
                     if remaining_edge(ind, 1) == current_min_vertex_index || remaining_edge(ind, 2) == current_min_vertex_index
                         remaining_edge(ind,:) = [];
                         tmp_num_deges = tmp_num_deges - 1;
                     else
                         ind =ind + 1;
                     end
                 end
            else
                % remove vertex from tmp_remaining_vertex
                tmp_remaining_vertex(:, current_min_vertex_col_ind_in_tmp_RV) = [];
            end
        end
        
    end
    
    
    
    
    % reorder the set index
    unusedindex =1; 
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
    num_subgroups_temp = set_count;
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
   
    
    number_of_subgroups = set_count;
    partition_result = set_node;
end