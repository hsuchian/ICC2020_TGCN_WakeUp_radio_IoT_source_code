% simulation2 main function
% generate nodes, input 

% max_range   : range for generating random points
% num_devices : number of devices in this group (same netmask)

function [SuccessProb_prim_mod, TransmissionTime_prim_mod, EnergyConsumption_prim_mod, Times_waked_prim_mod, partition_result_prim_mod, middle_point_prim_mod ] = S2_mainfunction_prim_final(max_range, num_devices, UAVradius, fscale, pointX, pointY, Length_of_hash)

   hash_inner = zeros(3,Length_of_hash);
   global hash1
   hash_inner(1,:) = hash1;
   global hash15
   hash_inner(2,:) = hash15;
   global hash2
   hash_inner(3,:) = hash2;
   
    % state parameters
    success_table = zeros(1, num_devices);
    times_of_being_waked = zeros(1, num_devices);
    Transmission_Time = zeros(1, num_devices);
    Energy_Consumption = zeros(1, num_devices);
    
    
    %Alpha = [0.49, 0.49, 0.34 , 0.24];
    Tslot = 4*10^(-3);    % slot duration = 4 ms, (in seconds)
    bitrate = 250*10^3;   % bit rate = 250 kbps
    WuCL = 44;            % wake up call length = 44 bits
    decodingTime = 0.3812*10^(-3);
    
    DeepE = 0.3*10^(-6);   % A
    LightE = 1.9*10^(-6);  % A
    ActiveE = 4.1*10^(-3); % A   Set V = 3
    
    
    % partition the devices, the index of subgroups are 1,2,3,....
    [partition_result_prim_mod, num_subgroups_our, middle_point_prim_mod] = S2_partition_prim_final_new_initial(pointX, pointY, UAVradius);
    
    
    % calculate real subgroups (),   some other points might be covered and 
    % perform scenario1 on each subgroup
    accumulate_time = 0;
    
    for subgroup_index = 1:num_subgroups_our 
        real_sbg_cur = [];        % indexes of devices in the current subgroup
        mid_x = middle_point_prim_mod(1, subgroup_index);
        mid_y = middle_point_prim_mod(2, subgroup_index);
        
        % construct the current real subgroup
        for device_index = 1:num_devices
            % calculate the distance between middle point 
            % if ( UAVradius^2 >= ((mid_x-pointX(device_index))^2 +  (mid_y - pointY(device_index))^2) )
            if (  UAVradius >=  norm( [mid_x, mid_y]- [pointX(device_index), pointY(device_index)] )   )
                real_sbg_cur = [real_sbg_cur, device_index];
            end
        end
        
        if isempty(real_sbg_cur ) 
            continue;
        end
        
        % perform scenario1 on real_sbg_cur
        % the mapping from real to temporary index is the index in real_sbg_cur
        numSensor_sbg = length(real_sbg_cur);
        FrameL = round(numSensor_sbg *fscale);
        
        % alpha is the collision probability
        if numSensor_sbg > Length_of_hash          % the number of devices exceeds the record, use the maximum record
            fprintf('number of devices exceeds : %d > %d (max record)\n', numSensor_sbg, Length_of_hash);
            %pause
            alpha = 1 - hash_inner(fscale*2-1, Length_of_hash);
        else
            alpha = 1 - hash_inner(fscale*2-1, numSensor_sbg);
        end
        
        
        RFrameM = round(alpha*numSensor_sbg); % alpha*N, where alpha is the collision probability of hash function~to fscale 
        %RFrameM = round(alpha*numSensor_sbg*fscale); 
        % [Success_table_sub, TransmissionTime_sub, EnergyConsumption_sub] = SimWakeUpRadio_for_S2(numSensor_sbg, FrameL, RFrameM, s_ind);
        [Success_table_sub, TransmissionTime_sub, EnergyConsumption_sub] = SimWakeUpRadio_for_S2_fast(numSensor_sbg, FrameL, RFrameM, alpha);
               
        
        % interpret the result, update the result
        for tmpInd= 1:length(real_sbg_cur)
            real_index = real_sbg_cur(tmpInd);
            
            if success_table(real_index) == 0 && Success_table_sub(tmpInd) == 1   % update the transmission time according to TransmissionTime_sub
                Transmission_Time(real_index) = Transmission_Time(real_index) + TransmissionTime_sub(tmpInd);
                
            elseif   success_table(real_index) == 0 && Success_table_sub(tmpInd) == 0   % failed
                Transmission_Time(real_index) = Transmission_Time(real_index) + (FrameL + RFrameM)*Tslot +  44 * (decodingTime ) ;
                if (FrameL + RFrameM)*Tslot +  44 * (decodingTime ) ~= TransmissionTime_sub(tmpInd)
                    fprintf('Transmission time Error at ours\n');
                    fprintf('   %d\n', (FrameL + RFrameM)*Tslot);
                    fprintf('   %d\n', TransmissionTime_sub(tmpInd));
                end
            % elseif  success_table(real_index) == 1 && Success_table_sub(tmpInd) == 0
                % the transmission time is fixed after successd 
            % elseif  success_table(real_index) == 1 && Success_table_sub(tmpInd) == 1 
                % the transmission time is fixed 
            end
                
            success_table(real_index) = ( (success_table(real_index) + Success_table_sub(tmpInd)) > 0);
            times_of_being_waked(real_index) = times_of_being_waked(real_index) + 1;                  
            Energy_Consumption(real_index) = Energy_Consumption(real_index) + EnergyConsumption_sub(tmpInd);
        end
        
        % devices not in range, thus, all in deep sleep mode
        deep_sleep_node = setdiff(1:num_devices, real_sbg_cur);
        for tmpInd= 1:length(deep_sleep_node)
            real_index = deep_sleep_node(tmpInd);
            
            if success_table(real_index) == 0                % if the device is not success yet, update the transmission time, plus the WuC decodeing time 
                Transmission_Time(real_index) = Transmission_Time(real_index) + (FrameL + RFrameM)*Tslot + 44 * (decodingTime );
            end
            Energy_Consumption(real_index) = Energy_Consumption(real_index) + (FrameL + RFrameM)*Tslot*DeepE*3 + 3*DeepE*(  44*( decodingTime   ) );    % always update the energy consumption of every other devices 
        end
        
    end  % end of for loop 
    
    %SuccessProb = success_table;
    %TransmissionTime = Transmission_Time;
    %EnergyConsumption = Energy_Consumption;
    
    Times_waked_prim_mod = mean(times_of_being_waked);
    SuccessProb_prim_mod = mean(success_table);
    % SuccessProb =times_of_being_waked;
    TransmissionTime_prim_mod = mean(Transmission_Time);
    EnergyConsumption_prim_mod = mean(Energy_Consumption);
end




