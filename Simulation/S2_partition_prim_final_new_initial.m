% scenario 2 
% Partition algorithm

function [partition_result, number_of_subgroups, middle_point] = S2_partition_prim_final_new_initial(X, Y, UAVradius)
    
    int_max = 2147483647;
    n_of_nodes = length(X);
    
    
    % use a 2*n matrix to perform the set operation (the set state),      [  ID ;  set ID ]
    set_node = [  1:length(X) ;  -1*ones(1, length(X))  ]; 
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
    remaining_edge = sortrows(remaining_edge, 3);       % the remaining edge set
    % Distance = sortrows(Distance, 3);   % sort in increasing order
  
    % remove the edges that has length > 3 UAVradius since it might be useless in our partition
    for remove_ind = 1:length( remaining_edge(:,1) )
        if remaining_edge(remove_ind, 3) > 3*UAVradius
            break;
        end
    end
    tmp_edge_num = length(remaining_edge(:,1));
    remaining_edge(remove_ind:tmp_edge_num, :) = [];
    
    
    % -------------------------------------- new partition algorithm ------------------------------
                     % [ vertex index ; distance to the current subgroup ; in current subgroup or not  ]
    remaining_vertex = [ 1:n_of_nodes ; int_max*ones(1, n_of_nodes)     ; zeros(1, n_of_nodes)  ];   
    
    %                  [ ori group index ; cen x ; cen y  ]
    middle_point_ori = [-1*ones(1, n_of_nodes)  ; -1*ones(1, n_of_nodes)  ];
    
    %fprintf('Start partitioning new mod prim\n');
    
    set_count = 0;
    while isempty(remaining_vertex)  ~= 1
        set_count = set_count+1;
        
        %[closet_vertex_last_round_val, closet_vertex_last_round_col_ind_in_RV] = min(remaining_vertex(2,:));
        
        %initialize remaining_vertex
        remaining_vertex(2,:) = int_max*ones(1, length(remaining_vertex(1,:)) );
        remaining_vertex(3,:) = zeros(1, length(remaining_vertex(1,:)) );
        

        % choose extremal point as initial 
        max_norm = 0;
        for ind_5 = 1:length(remaining_vertex(1,:))
            tmp_vertex_ind = remaining_vertex(1, ind_5);
            if norm([X(tmp_vertex_ind) Y(tmp_vertex_ind)]  - [100 100] ) > max_norm
                max_norm = norm([X(tmp_vertex_ind) Y(tmp_vertex_ind)]  - [100 100] );
                min_vertex_ind = tmp_vertex_ind;
            end
        end
            

        %min_vertex_ind = remaining_vertex(1, closet_vertex_last_round_col_ind_in_RV);
  
        
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
        remove_edges_row_ind_in_RE_out = find(remaining_edge(:,1) == min_vertex_ind | remaining_edge(:, 2) == min_vertex_ind);
        remaining_edge(remove_edges_row_ind_in_RE_out, :) = [];

        
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
                 remove_edges_row_ind_in_RE = find(remaining_edge(:,1) == current_min_vertex_index | remaining_edge(:, 2) == current_min_vertex_index);
                 remaining_edge(remove_edges_row_ind_in_RE, :) = [];
             
            else
                % remove vertex from tmp_remaining_vertex
                tmp_remaining_vertex(:, current_min_vertex_col_ind_in_tmp_RV) = [];
            end
        end  % end of inner while "while isempty(tmp_remaining_vertex) ~= 1"
        
        % find center of current set
        ind_of_current_subgroup = find(set_node(2,:) == set_count);
        [cen_x, cen_y] = find_center( X(ind_of_current_subgroup), Y(ind_of_current_subgroup), UAVradius);
        
        middle_point_ori(1, set_count) = cen_x;
        middle_point_ori(2, set_count) = cen_y;
        
        
        % ---------------------------------------- include the extra node in the current subgroup ------------------------------------------
        t_ind = 1;
        
        while t_ind <= length(remaining_vertex(1,:))
            
            tmp_ver_ind = remaining_vertex(1, t_ind);
            
            if norm([X(tmp_ver_ind), Y(tmp_ver_ind)] - [cen_x, cen_y]) <= UAVradius
                remaining_vertex(:, t_ind) = [];
                
                set_node(2, tmp_ver_ind) = set_count;
                set_element(set_count, n_of_nodes+1) = set_element(set_count, n_of_nodes+1) + 1;
                n_nodes_in_set_count = set_element(set_count, n_of_nodes+1);
                set_element(set_count, n_nodes_in_set_count) =  tmp_ver_ind;
                
                % remove the adjacent edges of this vertex
                remove_edges_row_ind_in_RE_last = find(remaining_edge(:,1) == tmp_ver_ind | remaining_edge(:, 2) == tmp_ver_ind);
                remaining_edge(remove_edges_row_ind_in_RE_last, :) = [];
                
            else 
                % calculate the distance of remaining vertex and current set center 
                remaining_vertex(2, t_ind) = norm([X(tmp_ver_ind), Y(tmp_ver_ind)] - [cen_x, cen_y]);
                t_ind = t_ind +1;
            end
        end
        
        % calculate the distance of remaining vertex and current set center 
        %{
        for t_ind = 1:length(remaining_vertex(1,:))
            tmp_ind = remaining_vertex(1, t_ind);
            if (norm(  [ X(tmp_ind)  Y(tmp_ind)] - [cen_x  cen_y] ) ~= remaining_vertex(2, t_ind))
                fprintf('ERROR, wrong dis calculation'); 
            end
            remaining_vertex(2, t_ind) = norm(  [ X(tmp_ind)  Y(tmp_ind)] - [cen_x  cen_y] );
        end
        %}
        
    end
    
   
    
    % reorder the set index
    unusedindex =1; 
    hash_table = -1*ones(1, 10000);   % -1 means no hit history, o.w., gives an unusedindex
    
    middle_point = -1000*ones(2, set_count);
    
    for iind = 1:length(X)
        cur_subg_index = set_node(2, iind);
        if hash_table( cur_subg_index ) == -1
            
            % fill center to final returned array
            middle_point(1, unusedindex) = middle_point_ori(1, cur_subg_index);
            middle_point(2, unusedindex) = middle_point_ori(2, cur_subg_index);
            
            hash_table(cur_subg_index) = unusedindex;
            set_node(2, iind) = unusedindex;
            unusedindex = unusedindex + 1;
        else
            set_node(2, iind) = hash_table(cur_subg_index);
        end
    end
    
    
    % --------------------------------------------------- calculate the middle point of each subgroup for tree partition -----------------------------------
    % index stard from 1
    %{
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
        
        
        [mid_x, mid_y] = find_center(X(index_of_devices_sbg), Y(index_of_devices_sbg), UAVradius);
        middle_point(:, subgroup_index) = [mid_x; mid_y];
    end
    %}
    
    number_of_subgroups = set_count;
    partition_result = set_node;
end