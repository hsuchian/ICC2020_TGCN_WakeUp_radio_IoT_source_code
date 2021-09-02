% Simulation controller
Sim_time = 10;
numSensor = 100;

limitScale = 4;

Success_prob = zeros(1, Sim_time); 
Success_prob_R = zeros(1, Sim_time); 

scaleArray = 1:0.5:limitScale;

% Final result for each "scale"
Success= zeros(1, (limitScale-1)/0.5);
Success_R= zeros(1, (limitScale-1)/0.5);


index = 1;
for scale = scaleArray 
    
    for i = 1:Sim_time
        % SimWakeUpRadio(numSensor, FrameL, RFrameL)
        Success_prob(i) = SimHash(numSensor, round(numSensor*scale), round(numSensor/5));
        Success_prob_R(i) = SimRandom(numSensor, round(numSensor*scale), round(numSensor/5));
    end
    Success(index) = mean(Success_prob);
    Success_R(index) = mean(Success_prob_R);
    index = index + 1;
    Success_prob = zeros(1, Sim_time); 
    Success_prob_R = zeros(1, Sim_time); 
end


save data\SuccessProb