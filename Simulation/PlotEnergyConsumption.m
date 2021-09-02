% Plot 
clear;
clc;
close all;

Sim_time = 150;
num_devices = 50;
load(['data/S1_simtime_', num2str(Sim_time), '_numD_', num2str(num_devices)]);

f1 = figure(1);

hold on
decodingTime = 0.3812*10^(-3);

DeepE = 0.3*10^(-6);           % A
LightE = 1.9*10^(-6);          % A
ActiveE = 4.1*10^(-3);         % A   Set V = 3



alpha1 = 1 - hash1;
alpha1 = alpha1(round(scaleArray*numSensor));

alpha15 = 1 - hash15;
alpha15 = alpha15(round(scaleArray*numSensor));

alpha2 = 1 - hash2;
alpha2 = alpha2(round(scaleArray*numSensor));



Y = zeros(1, 9);
P1_1 = zeros(1, 9);
P1_2 = zeros(1, 9);
P1_3 = zeros(1, 9);

for i = 1:9
    P1_1(i) = 1 - alpha1(i);
    P1_2(i) = alpha1(i) * (1 - ( 1 / ( alpha1(i)*scaleArray(i)*numSensor ) ) ) ^ ( alpha1(i)*scaleArray(i)*numSensor - 1); 
    P1_3(i) = -P1_2(i) + alpha1(i); 
    
    Y(i) = P1_1(i) * ( ActiveE*3  +  3*LightE*(scaleArray(i)*numSensor-1)/2  ...
                                  +  3*DeepE*(scaleArray(i)*numSensor-1)/2   ...
                                  +  3*DeepE*alpha1(i)*scaleArray(i)*numSensor) ...
         + P1_2(i) * (2*ActiveE*3 +  3*LightE*(scaleArray(i)*numSensor-1) ...
                                  +  3*LightE*(alpha1(i)*scaleArray(i)*numSensor-1)/2 ...
                                  +  3*DeepE*(alpha1(i)*scaleArray(i)*numSensor-1)/2) ...
         + P1_3(i) * (2*ActiveE*3 +  3*LightE*(scaleArray(i)*numSensor-1) ...
                                  +  3*LightE*(alpha1(i)*scaleArray(i)*numSensor-1)/2 ...
                                  +  3*DeepE*(alpha1(i)*scaleArray(i)*numSensor-1)/2);
    
    %Y(i) = (1-alpha(1))*(ActiveE*3 + 3*LightE*(scaleArray(i)*numSensor-1))/2  +  alpha(1)*(LightE*3*(alpha(1)*scaleArray(i)*numSensor + 2*scaleArray(i)*numSensor-3) + 4*ActiveE*3)/2;
end

Y15 = zeros(1, 9);
P15_1 = zeros(1, 9);
P15_2 = zeros(1, 9);
P15_3 = zeros(1, 9);
for i = 1:9
    P15_1(i) = 1 - alpha15(i);
    P15_2(i) = alpha15(i) * (1 - ( 1 / (alpha15(i)*scaleArray(i)*numSensor) ) ) ^ ( alpha15(i)*scaleArray(i)*numSensor - 1); 
    P15_3(i) = -P15_2(i) + alpha15(i); 
    
    Y15(i) = P15_1(i) * ( ActiveE*3  +  3*LightE*(1.5*scaleArray(i)*numSensor-1)/2  ...
                                     +  3*DeepE*(1.5*scaleArray(i)*numSensor-1)/2   ...
                                     +  3*DeepE*(alpha15(i)*scaleArray(i)*numSensor)) ...
         + alpha15(i)* (2*ActiveE*3  +  3*LightE*(1.5*scaleArray(i)*numSensor-1) ...
                                     +  3*LightE*(alpha15(i)*scaleArray(i)*numSensor-1)/2 ...
                                     +  3*DeepE*(alpha15(i)*scaleArray(i)*numSensor-1)/2);    
    
    %{
    Y15(i) = P15_1(i) * ( ActiveE*3  +  3*LightE*(1.5*scaleArray(i)*numSensor-1)/2  ...
                                     +  3*DeepE*(1.5*scaleArray(i)*numSensor-1)/2   ...
                                     +  3*DeepE*(alpha15(i)*scaleArray(i)*numSensor)) ...
         + P15_2(i) * (2*ActiveE*3   +  3*LightE*(1.5*scaleArray(i)*numSensor-1) ...
                                     +  3*LightE*(alpha15(i)*scaleArray(i)*numSensor-1)/2 ...
                                     +  3*DeepE*(alpha15(i)*scaleArray(i)*numSensor-1)/2) ...
         + P15_3(i) * (2*ActiveE*3   +  3*LightE*(1.5*scaleArray(i)*numSensor-1) ...
                                     +  3*LightE*(alpha15(i)*scaleArray(i)*numSensor-1)/2 ...
                                     +  3*DeepE*(alpha15(i)*scaleArray(i)*numSensor-1)/2  );
%}
    %Y15(i) = (1-alpha(2))*(ActiveE*3 + 3*LightE*(1.5*scaleArray(i)*numSensor-1))/2  +  alpha(2)*(LightE*3*(alpha(2)*1.5*scaleArray(i)*numSensor + 2*1.5*scaleArray(i)*numSensor-3) + 4*ActiveE*3)/2;
end


Y2 = zeros(1, 9);
P2_1 = zeros(1, 9);
P2_2 = zeros(1, 9);
P2_3 = zeros(1, 9);
for i = 1:9
    P2_1(i) = 1 - alpha2(i);
    P2_2(i) = alpha2(i) * (1 - ( 1 / ( alpha2(i)*scaleArray(i)*numSensor ) ) ) ^ ( alpha2(i)*scaleArray(i)*numSensor - 1); 
    P2_3(i) = -P2_2(i) + alpha2(i); 
    
    Y2(i) = P2_1(i) * ( ActiveE*3  +  3*LightE*(2*scaleArray(i)*numSensor-1)/2  ...
                                  +  3*DeepE*(2*scaleArray(i)*numSensor-1)/2   ...
                                  +  3*DeepE*alpha2(i)*scaleArray(i)*numSensor) ...
         + alpha2(i) * (2*ActiveE*3 +  3*LightE*(2*scaleArray(i)*numSensor-1) ...
                                  +  3*LightE*(alpha2(i)*scaleArray(i)*numSensor-1)/2 ...
                                  +  3*DeepE*(alpha2(i)*scaleArray(i)*numSensor-1)/2);
    
    %{
    Y2(i) = P2_1(i) * ( ActiveE*3  +  3*LightE*(2*scaleArray(i)*numSensor-1)/2  ...
                                  +  3*DeepE*(2*scaleArray(i)*numSensor-1)/2   ...
                                  +  3*DeepE*alpha2(i)*scaleArray(i)*numSensor) ...
         + P2_2(i) * (2*ActiveE*3 +  3*LightE*(2*scaleArray(i)*numSensor-1) ...
                                  +  3*LightE*(alpha2(i)*scaleArray(i)*numSensor-1)/2 ...
                                  +  3*DeepE*(alpha2(i)*scaleArray(i)*numSensor-1)/2) ...
         + P2_3(i) * (2*ActiveE*3 +  3*LightE*(2*scaleArray(i)*numSensor-1) ...
                                  +  3*LightE*(alpha2(i)*scaleArray(i)*numSensor-1)/2 ...
                                  +  3*DeepE*(alpha2(i)*scaleArray(i)*numSensor-1)/2);
    
    %}
    %Y2(i) = (1-alpha(3))*(ActiveE*3 + 3*LightE*(2*scaleArray(i)*numSensor-1))/2  +  alpha(3)*(LightE*3*(alpha(3)*2*scaleArray(i)*numSensor + 2*2*scaleArray(i)*numSensor-3) + 4*ActiveE*3)/2;
end



Y = Y*(4*10^(-3));
Y = Y+  44 * (decodingTime ) * LightE*3 ;

Y15 = Y15*(4*10^(-3));
Y15 = Y15+  44 * (decodingTime ) * LightE*3 ;

Y2 = Y2*(4*10^(-3));
Y2 = Y2+  44 * (decodingTime ) * LightE*3 ;

%{
plot( scaleArray*numSensor , Energy_Consumption(1,:), 'b*');
plot( scaleArray*numSensor ,Y , 'b-.');
plot( scaleArray*numSensor , Energy_ConsumptionR(1,:), 'b-*');

plot( scaleArray*numSensor , Energy_Consumption(2,:), 'gs');
plot( scaleArray*numSensor ,Y15 , 'g--');
plot( scaleArray*numSensor , Energy_ConsumptionR(2,:), 'g-s');

plot( scaleArray*numSensor , Energy_Consumption(3,:), 'rd');
plot( scaleArray*numSensor ,Y2 , 'r:');
plot( scaleArray*numSensor , Energy_ConsumptionR(3,:), 'r-d');
%}

plot( scaleArray*numSensor , Energy_Consumption(1,:), 'bh');
plot( scaleArray*numSensor ,Y , 'b:', 'LineWidth', 1.2);
plot( scaleArray*numSensor , Energy_ConsumptionR(1,:), 'b-o');

plot( scaleArray*numSensor , Energy_Consumption(2,:), 'gs');
plot( scaleArray*numSensor ,Y15 , 'g--');
plot( scaleArray*numSensor , Energy_ConsumptionR(2,:), 'g-d');

plot( scaleArray*numSensor , Energy_Consumption(3,:), 'r^');
plot( scaleArray*numSensor ,Y2 , 'r-.');
plot( scaleArray*numSensor , Energy_ConsumptionR(3,:), 'r-v');

%Y = (1- 1./scaleArray/numSensor).^(numSensor-1);
%plot( scaleArray , Y, 'b--');
%plot( scaleArray , Success, 'g:^');

%lgd = legend('Ours L = N', 'Ours (ana) L = N \alpha = 0.49',  'CSMA L = N', 'Ours L = 1.5N', 'Ours (ana) L = 1.5N \alpha = 0.34', 'CSMA L = 1.5N', 'Ours L = 2N', 'Ours (ana) L = 2N, \alpha = 0.24', 'CSMA L =2N');
lgd = legend('Ours L = N', 'Ours L = N (ana)',  'CSMA L = N', 'Ours L = 1.5N', 'Ours L = 1.5N (ana)', 'CSMA L = 1.5N', 'Ours L = 2N', 'Ours L = 2N (ana)', 'CSMA L = 2N');
lgd.Position(1) = 0.17;
lgd.Position(2) = 0.55;
lgd.FontSize = 9;

xticks(100:10:200);

set(gca,'FontSize',12);

% lgd = legend('Ours L = N','Ours L = 1.5N','Ours L =2N', 'CSMA L = N','CSMA L = 1.5N','CSMA L =2N');
%set(lgd,'Interpreter','latex');
%lgd.FontSize=12;
%lgd.FontWeight='bold';


xlabel('Number of devices');  
ylabel('Energy consumption (Joule)');
%ylim([0.5*10^(-4) 1*10^(-4)]);
hold off;
%----------------------------------------- sub plot -------------------------------


axes('Position',[0.65 0.27 0.2 0.2]);

hold on;

plot( scaleArray*numSensor , Energy_Consumption(1,:), 'bh');
plot( scaleArray*numSensor ,Y , 'b:', 'LineWidth', 1.2);

plot( scaleArray*numSensor , Energy_Consumption(2,:), 'gs');
plot( scaleArray*numSensor ,Y15 , 'g--');

plot( scaleArray*numSensor , Energy_Consumption(3,:), 'r^');
plot( scaleArray*numSensor ,Y2 , 'r-.');


annotation('textarrow',[0.595 0.635], [0.14 0.235]);
%annotation('arrow', [0 0], [4 4]);
xlim([160 200]);
ylim([6*10^(-5) 8.5*10^(-5)]);
hold off;

%figure1 = figure(1);
%set(figure1, 'PaperUnits', 'inches');
%set(figure1, 'PaperPosition', [0 0 4.7 5]);

% set position
%set(figure1, 'Units', 'inches');
%set(figure1, 'Position', [0 0 8 8]);

print(f1, '-depsc', ['picture\','S1_EnergyConsumption_Sim_', num2str(Sim_time),'_N_', num2str(num_devices),  '.eps']);
print(f1, '-dpng', ['picture\','S1_EnergyConsumption_Sim_', num2str(Sim_time),'_N_', num2str(num_devices),  '.png']);