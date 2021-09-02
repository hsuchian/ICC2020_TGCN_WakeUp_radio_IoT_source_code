
% wake up radio simulation 
function  [SuccessProb, TransmissionTime, EnergyConsumption] = SimWakeUpRadio(numSensor, FrameL, RFrameM, maxseed)

    % sim_time = 1;

    % Scenario 1
    % Area_length_X  = 100;  % m
    % Area_length_Y = 100;   % m
    % numSensor = 10;
    % FrameL = 1.5*numSensor;
    % RFrameL = 20;

    % Generate the location of sensors
    % Xaxis = Area_length_X * rand(1,numSensor);
    % Yaxis = Area_length_X * rand(1,numSensor);

    % generate the ID of each sensor
    ID_int = 1:numSensor;
    ID_bi = de2bi(ID_int, 8);
    threshold = 0.9;
    seed_rand_index = randi([-100000, 100000]);
    limit_iteration = 200;
    
    Energy_consumption = zeros(1, numSensor);    % (ID) [ deep,  light, active   ] 
    Transmission_Time  = zeros(1,  numSensor );    % (ID) [ deep time ,  light time , active time ]   in time slot, 
    
    % Parameters
    Tslot = 4*10^(-3);    % slot duration = 4 ms, (in seconds)
    bitrate = 250*10^3;   % bit rate = 250 kbps
    WuCL = 44;            % wake up call length = 44 bits
    decodingTime = 0.3812*10^(-3);
    
    DeepE = 0.3*10^(-6);   % A
    LightE = 1.9*10^(-6);  % A
    ActiveE = 4.1*10^(-3); % A   Set V = 3

    % Simulation start    
    State_of_sensor = zeros(1, numSensor);            % 0:failed, 1:success
    maxValue = 0;
    max_state_of_sensor = zeros(1, numSensor);            % 0:failed, 1:success
    
    iter = 0;
    % while maxValue < threshold
    if maxseed ~= 0 
        hashResult = zeros(1, numSensor);                 % hash table
        HashTotalResult = zeros(1,FrameL);
        State_of_sensor = zeros(1, numSensor);            % 0:failed, 1:success
        %seed = randi([0 500]);
        seed = maxseed;
        for i = 1:numSensor
            hashResult(i) = Hash(ID_int(i), FrameL, seed);
            % fprintf('%d \n', hashResult(i));
            HashTotalResult( hashResult(i) ) =  HashTotalResult(hashResult(i) ) + 1;
        end
    
        for i = 1:numSensor
            if HashTotalResult( hashResult(i) ) == 1
                State_of_sensor(i) = 1; 
            end
        end
        if mean( State_of_sensor) > maxValue
            maxValue = mean( State_of_sensor);
            max_state_of_sensor = State_of_sensor;
            max_hash_result = hashResult;         % store the hash table
        end
        iter = iter +1;
        if iter > limit_iteration
            iter = 1;
        end
    end
    
    % fprintf('%d', maxValue);
    
    % hash table : max_hash_result
    % calculate the energy consumption and transmission time for success users (based on hash table : max_hash_result)
    for i = 1:numSensor
        if max_state_of_sensor(i) == 1
            Transmission_Time(i) = max_hash_result(i)*Tslot;
            Energy_consumption(i) = 3*DeepE*(FrameL - max_hash_result(i)+ RFrameM)*Tslot +  3*LightE*(max_hash_result(i) - 1)*Tslot +   3*ActiveE*Tslot;               % (ID) [ deep,  light, active   ] 
        end 
    end
    
    
    
    % unsuccess sensor : max_state_of_sensor
    if RFrameM == 0
        SuccessProb = maxValue; 
        TransmissionTime = 0;
        EnergyConsumption = 0;
        return;
    end
    % ----------------------------------------------- Performa random frame---------------------------------------------------------------
    random_choice = zeros(1,numSensor);
    RframeResult = zeros(1, RFrameM);
    State_of_sensor_final = max_state_of_sensor;
    for i = 1:numSensor
        if max_state_of_sensor(i) == 0
           random_choice(i) = randi([1 RFrameM]);
           RframeResult(random_choice(i)) = RframeResult(random_choice(i)) + 1;
        end
    end
    
    for i = 1:numSensor
        if max_state_of_sensor(i) == 0
            if  RframeResult(random_choice(i)) == 1
                State_of_sensor_final(i) = 1;
            end
        end
    end
    SuccessProb = mean(State_of_sensor_final);
    
    
    % calculate the energy consumption of failed-1 devices
    for i = 1:numSensor
        if max_state_of_sensor(i) == 0 && State_of_sensor_final(i) == 1
            if random_choice(i) == 0
                fprintf('ERROR in random result\n');
            end
            Transmission_Time(i) =  (FrameL + random_choice(i) )*Tslot;              % in terms of time slot
            Energy_consumption(i) = 3*DeepE*(RFrameM - random_choice(i) )*Tslot  +   3*LightE*(FrameL -1 + random_choice(i) -1)*Tslot +   3*2*ActiveE*Tslot;               % (ID) [ deep,  light, active   ] 
        end 
    end
    
    
    for i = 1:numSensor
        if State_of_sensor_final(i) == 0
            
            Transmission_Time(i) =  (FrameL + RFrameM)*Tslot;              % in terms of time slot
            Energy_consumption(i) = 3*DeepE*(RFrameM - random_choice(i) )*Tslot  +   3*LightE*(FrameL -1 + random_choice(i) -1)*Tslot +   3*2*ActiveE*Tslot;               % (ID) [ deep,  light, active   ] 
        end 
    end
    
    TransmissionTime = mean(Transmission_Time)   +  44 * (decodingTime ) ; % plus wake up call decoding time 44bits * 0.3812 ms
    EnergyConsumption = mean(Energy_consumption)  +  3*LightE*(  44*( decodingTime   ) );
    
end