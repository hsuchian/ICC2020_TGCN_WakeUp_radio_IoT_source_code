function  [SuccessProb, maxseed] = SimHash(numSensor, FrameL, RFrameL)
    % wake up radio simulation 

    % sim_time  1;

    % Scenario 1
    Area_length_X  = 100; % m
    Area_length_Y = 100;   % m
    threshold = 0.9;
    limit_iteration = 1000;
    % numSensor = 10;
    % FrameL = 1.5*numSensor;
    % RFrameL = 20;

    % Generate the location of sensors
    Xaxis = Area_length_X * rand(1,numSensor);
    Yaxis = Area_length_X * rand(1,numSensor);

    % generate the ID of each sensor
    ID_int = 1:numSensor;
    ID_bi = de2bi(ID_int, 8);

    % Result data

    State_of_sensor = zeros(1, numSensor);            % 0:failed, 1:success
    maxValue = 0;
    iter = 0;
    % Simulation start
    while maxValue < threshold
        hashResult = zeros(1, numSensor);                 % hash table
        HashTotalResult = zeros(1,FrameL);
        State_of_sensor = zeros(1, numSensor);            % 0:failed, 1:success
        %seed = randi([0 500]);
        seed = iter - 1000;
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
            maxseed = seed;
        end
        iter = iter +1;
        if iter > limit_iteration
            break;
        end
    end
    % ----------------------------------------------- Performa random frame---------------------------------------------------------------
 %{ 
    random_choice = zeros(1,numSensor);
    RframeResult = zeros(1, RFrameL);
    State_of_sensor_final = State_of_sensor;
    for i = 1:numSensor
        if State_of_sensor(i) == 0
           random_choice(i) = randi([1 RFrameL]);
           RframeResult(random_choice(i)) = RframeResult(random_choice(i)) + 1;
        end
    end
    
    for i = 1:numSensor
        if State_of_sensor_final(i) == 0
            if  RframeResult(random_choice(i)) == 1
                State_of_sensor_final(i) = 1;
            end
        end
    end
    %}
    SuccessProb = maxValue;
end