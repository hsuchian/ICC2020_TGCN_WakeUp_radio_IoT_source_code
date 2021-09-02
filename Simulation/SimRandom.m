function  SuccessProb = SimRandom(numSensor, FrameL, RFrameL)
    % wake up radio simulation 

    % sim_time = 1;

    % Scenario 1
    Area_length_X  = 100; % m
    Area_length_Y = 100;   % m
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


    % Simulation start
    hashResult = zeros(1, numSensor);                 % hash table
    HashTotalResult = zeros(1,FrameL);
    State_of_sensor = zeros(1, numSensor);            % 0:failed, 1:success
    seed = randi([0 500]);
    
    
    for i = 1:numSensor
        % hashResult(i) = Hash(ID_int(i), FrameL, seed);
        hashResult(i) = randi([1 FrameL]);
        % fprintf('%d \n', hashResult(i));
        HashTotalResult( hashResult(i) ) =  HashTotalResult(hashResult(i) ) + 1;
    end
    
    for i = 1:numSensor
        if HashTotalResult( hashResult(i) ) == 1
             State_of_sensor(i) = 1; 
        end
    end
    
    SuccessProb = mean(State_of_sensor);
end