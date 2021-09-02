function  [SuccessProb, record, rrr] = SimHash1(numSensor, FrameL, seed_start, average_num)

    threshold = 0.9;
    limit_iteration = 100000;

    % generate the ID of each sensor
    ID_int = 1:numSensor;

    % Result data

    State_of_sensor = zeros(1, numSensor);            % 0:failed, 1:succes
    
    record = zeros(2, limit_iteration);
    
    seed = seed_start;
    maxValue = 0;
    
    
    % Simulation start
    for iter = 1:limit_iteration 
        hashResult = zeros(1, numSensor);                 % hash table
        HashTotalResult = zeros(1,FrameL);
        State_of_sensor = zeros(1, numSensor);            % 0:failed, 1:success
        seed = seed_start + iter;
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
        
        record(1, iter) = seed_start + iter;
        record(2, iter) = mean( State_of_sensor);
       
    end
    record = (sortrows(record.', -2)).';
    SuccessProb = maxValue;
    rrr = mean(record(2,1:average_num));
    %fprintf('x = %d, %d\n',numSensor ,mean(record(2,1:10)));
end