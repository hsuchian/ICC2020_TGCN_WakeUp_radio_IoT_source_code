% simulation2 main function
% generate nodes, input 

% max_range   : range for generating random points
% num_devices : number of devices in this group (same netmask)

function [SuccessProb_sq, TransmissionTime_sq, EnergyConsumption_sq, Times_waked_sq, partition_result_sq, middle_point_sq] = S2_mainfunction_sq(max_range, num_devices, UAVradius, fscale, pointX, pointY, Length_of_hash)

   hash_inner = zeros(3,Length_of_hash);
   global hash1
   hash_inner(1,:) = hash1;
   global hash15
   hash_inner(2,:) = hash15;
   global hash2
   hash_inner(3,:) = hash2;
   
    % state parameters
    success_table_sq = zeros(1, num_devices);
    times_of_being_waked_sq = zeros(1, num_devices);
    Transmission_Time_sq = zeros(1, num_devices);
    Energy_Consumption_sq = zeros(1, num_devices);
    
    
    % Alpha = [0.49, 0.49, 0.34 , 0.24];
    Tslot = 4*10^(-3);    % slot duration = 4 ms, (in seconds)
    bitrate = 250*10^3;   % bit rate = 250 kbps
    WuCL = 44;            % wake up call length = 44 bits
    decodingTime = 0.3812*10^(-3);
    
    DeepE = 0.3*10^(-6);   % A
    LightE = 1.9*10^(-6);  % A
    ActiveE = 4.1*10^(-3); % A   Set V = 3
   
    
    % ===================================================== start trivial partition method (square) ======================================================
    
    
     [partition_result_sq, num_subgroups_square, midpoint_sq_1 ] = S2_EvenPartition_square(pointX, pointY, UAVradius, max_range);
    
    % calculate the middle point of each subgroup for tree partition
    middle_point_sq = midpoint_sq_1;            

                                                             
    if num_subgroups_square > (ceil(max_range/( ( 2^(1/2) )*UAVradius  ) ))^2
        fprintf('Error for calculating subgroup boundaries and the number\n');
        pause  
    end
    
    % calculate real subgroups (),   some other points might be covered and 
    % perform scenario1 on each subgroup
    for subgroup_index = 1:num_subgroups_square
        real_sbg_cur = [];        % indexes of devices in the current subgroup
        mid_x = middle_point_sq(1, subgroup_index);
        mid_y = middle_point_sq(2, subgroup_index);
        
        % construct the current real subgroup
        for device_index = 1:num_devices
            % calculate the distance between middle point 
            if (  UAVradius >=  norm( [mid_x, mid_y]- [pointX(device_index), pointY(device_index)] )   )
            %if ( UAVradius^2 >= (  (mid_x- pointX(device_index))^2  +  (mid_y - pointY(device_index))^2)  )
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
            % pause
            alpha = 1 - hash_inner(fscale*2-1, Length_of_hash);
        else
            alpha = 1 - hash_inner(fscale*2-1, numSensor_sbg);
        end
        
        RFrameM = round(alpha*numSensor_sbg); % alpha*N, where alpha is the collision probability of hash function~to fscale 
        %if RFrameM == 0
        %    fprintf('');
        
        % [Success_table_sq_sub, TransmissionTime_sq_sub, EnergyConsumption_sq_sub] = SimWakeUpRadio_for_S2(numSensor_sbg, FrameL, RFrameM, s_ind);
        [Success_table_sq_sub, TransmissionTime_sq_sub, EnergyConsumption_sq_sub] = SimWakeUpRadio_for_S2_fast(numSensor_sbg, FrameL, RFrameM, alpha);
        
        % interpret the result, update the result
        for tmpInd= 1:length(real_sbg_cur)
            real_index = real_sbg_cur(tmpInd);
            
            if success_table_sq(real_index) == 0 && Success_table_sq_sub(tmpInd) == 1   % update the transmission time according to TransmissionTime_sub
                Transmission_Time_sq(real_index) = Transmission_Time_sq(real_index) + TransmissionTime_sq_sub(tmpInd);
                
            elseif   success_table_sq(real_index) == 0 && Success_table_sq_sub(tmpInd) == 0   % failed
                Transmission_Time_sq(real_index) = Transmission_Time_sq(real_index) + (FrameL + RFrameM)*Tslot +  44 * (decodingTime ) ;
                if (FrameL + RFrameM)*Tslot +  44 * (decodingTime ) ~= TransmissionTime_sq_sub(tmpInd)
                    fprintf('Transmission time Error at sq\n');
                    fprintf('   %d\n', (FrameL + RFrameM)*Tslot);
                    fprintf('   %d\n', TransmissionTime_sub(tmpInd));
                end
            % elseif  success_table_sq(real_index) == 1 && Success_table_sq_sub(tmpInd) == 0
                % the transmission time is fixed 
            % elseif  success_table_sq(real_index) == 1 && Success_table_sq_sub(tmpInd) == 1 
                % the transmission time is fixed 
            end
                        
            success_table_sq(real_index) = ( (success_table_sq(real_index) + Success_table_sq_sub(tmpInd)) > 0);
            
            % change this "tims_of_being_waked" into transmission time, also recalculate the energy consumption
            times_of_being_waked_sq(real_index) = times_of_being_waked_sq(real_index) + 1;
            Energy_Consumption_sq(real_index) = Energy_Consumption_sq(real_index) + EnergyConsumption_sq_sub(tmpInd);
        end
        
        % devices not in the range, thus, all of them are in deep sleep mode
        deep_sleep_node = setdiff(1:num_devices, real_sbg_cur);
        for tmpInd= 1:length(deep_sleep_node)
            real_index = deep_sleep_node(tmpInd);
            if success_table_sq(real_index) == 0                % if the device is not success yet, update the transmission time, plus decoding time 
                Transmission_Time_sq(real_index) = Transmission_Time_sq(real_index) + (FrameL + RFrameM)*Tslot + 44 * (decodingTime );
            end
            Energy_Consumption_sq(real_index) = Energy_Consumption_sq(real_index) + (FrameL + RFrameM)*Tslot*DeepE*3 + 3*DeepE*(  44*( decodingTime   ) );
        end
        
    end  % end of for loop 
    
    %SuccessProb_sq = success_table_sq;
    %TransmissionTime_sq = Transmission_Time_sq;
    %EnergyConsumption_sq = Energy_Consumption_sq;
    Times_waked_sq = mean(times_of_being_waked_sq);
    SuccessProb_sq = mean(success_table_sq);
    % SuccessProb_sq =times_of_being_waked_sq;
    TransmissionTime_sq = mean(Transmission_Time_sq);
    EnergyConsumption_sq = mean(Energy_Consumption_sq);
    
end




