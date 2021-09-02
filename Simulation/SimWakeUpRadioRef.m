
% simulation for CSMA based contention


function  [SuccessProb, TransmissionTime, EnergyConsumption] = SimWakeUpRadioRef(numSensor, FrameL, RFrameM)
% function  [SuccessProb] = SimWakeUpRadioRef(numSensor, FrameL, RFrameM)    
    % Parameters
    Tslot = 4*10^(-3);     % slot duration = 4 ms, (in seconds)
    bitrate = 250*10^3;    % bit rate = 250 kbps
    WuCL = 20;             % wake up call length = 20 bits
    
    DeepE = 0.3*10^(-6);           % A
    LightE = 1.9*10^(-6);          % A
    ActiveE = 4.1*10^(-3);         % A   Set V = 3
    decodingTime = 0.3812*10^(-3);
    
    CWsize = 32;
    CWslot = 320* 10^(-6); % 320 micro seconds
    
    ID = 1:numSensor;
    State_of_sensor = zeros(1, numSensor);            % 0:failed, 1:success
    Transmission_time = zeros(1, numSensor);  
    Energy_consumption = zeros(1, numSensor); 
    collisiontime = zeros(1, numSensor); 
    
    for time = 1: FrameL + RFrameM     % each round
        Backoff = (CWsize+5)*ones(1, numSensor);
        ContentionResult = zeros(1, CWsize);         % final result
        
        for i = 1:numSensor 
            if State_of_sensor(i) == 0
                Backoff(i) = randi(CWsize);
                ContentionResult(Backoff(i)) = ContentionResult(Backoff(i)) +1 ;
            end
        end
        
        [minBackoff, minIndice] = min(Backoff);
        if minBackoff > CWsize
            break;
        end
        if ContentionResult(minBackoff) > 1
            x = find(Backoff == minBackoff);
            for j = 1:length(x)
                collisiontime(x(j)) = collisiontime(x(j)) +1;
            end
            continue;
        else  % success
            State_of_sensor(minIndice) = 1;
            Transmission_time(minIndice) = time*(Tslot + CWslot*CWsize ) + decodingTime*WuCL;
            Energy_consumption(i) = (Tslot + CWslot*CWsize )*(collisiontime(i)+1)*3*ActiveE   + ... 
                                    (Tslot + CWslot*CWsize )*(time - collisiontime(i) - 1)*3*LightE  + ...
                                    (Tslot + CWslot*CWsize )*(FrameL + RFrameM - time)*3*DeepE  + ...
                                    3*LightE*(  WuCL*( decodingTime   ) ) ;
        end
    end
    
    
    % calculate the failed sensor 
    for i = 1:numSensor 
        if State_of_sensor(i) == 0
            Transmission_time(i) = (FrameL + RFrameM)*(Tslot + CWslot*CWsize ) + decodingTime*WuCL;
            Energy_consumption(i) = (Tslot + CWslot*CWsize )*(collisiontime(i))*3*ActiveE   + ... 
                                    (Tslot + CWslot*CWsize )*(FrameL + RFrameM - collisiontime(i))*3*LightE  + ...
                                    3*LightE*(  WuCL*( decodingTime   ) ) ;
        end
    end
    
    TransmissionTime = mean(Transmission_time);
    EnergyConsumption = mean(Energy_consumption);
    SuccessProb = mean(State_of_sensor);
end