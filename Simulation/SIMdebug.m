% Simulation controller
Sim_time = 100;
numSensor = 50;



scaleArray = 1:0.25:4;
scaleArray_frame_length = 1:0.5:2;
Success= zeros(1, (2-1)/0.25 + 1);

Success_Prob = zeros(3, (4-1)/0.25 +1); 
Transmission_Time = zeros(3, (4-1)/0.25 +1); 
Energy_Consumption = zeros(3, (4-1)/0.25 +1); 

Success_ProbR = zeros(3, (4-1)/0.25 +1); 
Transmission_TimeR = zeros(3, (4-1)/0.25 +1); 
Energy_ConsumptionR = zeros(3, (4-1)/0.25 +1); 
alpha = [0.4, 0.24, 0.15];;

index1 = 1;

for fscale = scaleArray_frame_length
    index2 = 1;
    for scale = scaleArray 
    
        for i = 1:Sim_time
            % SimWakeUpRadio(numSensor, FrameL, RFrameL)
            [SuccessProbR, TransmissionTimeR, EnergyConsumptionR] = SimWakeUpRadioRef(round(numSensor*scale), round(fscale*numSensor*scale), round(alpha(index1)*numSensor*scale));
            % [SuccessProb, TransmissionTime, EnergyConsumption] = SimWakeUpRadio(round(numSensor*scale), round(fscale*numSensor*scale), round(alpha(index1)*numSensor*scale));
            
            SuccessProb=0;
            TransmissionTime=0;
            EnergyConsumption =0;
            Success_Prob(index1, index2) = Success_Prob(index1, index2) + SuccessProb;
            Transmission_Time(index1, index2) = Transmission_Time(index1, index2) + TransmissionTime;
            Energy_Consumption(index1, index2) = Energy_Consumption(index1, index2) + EnergyConsumption;
            
            Success_ProbR(index1, index2) = Success_ProbR(index1, index2) + SuccessProbR;
            Transmission_TimeR(index1, index2) = Transmission_TimeR(index1, index2) + TransmissionTimeR;
            Energy_ConsumptionR(index1, index2) = Energy_ConsumptionR(index1, index2) + EnergyConsumptionR;
        end
        Success_Prob(index1, index2) = Success_Prob(index1, index2)/Sim_time;
        Transmission_Time(index1, index2) = Transmission_Time(index1, index2)/Sim_time;
        Energy_Consumption(index1, index2) = Energy_Consumption(index1, index2)/Sim_time;
        
        Success_ProbR(index1, index2) = Success_ProbR(index1, index2)/Sim_time;
        Transmission_TimeR(index1, index2) = Transmission_TimeR(index1, index2)/Sim_time;
        Energy_ConsumptionR(index1, index2) = Energy_ConsumptionR(index1, index2)/Sim_time;
        index2 = index2 + 1;
    end
    index1 = index1 + 1;
end


save data\netscaletest