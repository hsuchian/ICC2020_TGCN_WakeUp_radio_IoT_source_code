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


alpha1 = 1 - hash1;
alpha1 = alpha1(round(scaleArray*numSensor));

alpha15 = 1 - hash15;
alpha15 = alpha15(round(scaleArray*numSensor));

alpha2 = 1 - hash2;
alpha2 = alpha2(round(scaleArray*numSensor));




% analysis
Y = zeros(1, 9);
P1_1 = zeros(1, 9);
P1_2 = zeros(1, 9);
P1_3 = zeros(1, 9);
for i = 1:9
    P1_1(i) = 1 - alpha1(i);
    P1_2(i) = alpha1(i) * (1 - ( 1 / ( alpha1(i)*scaleArray(i)*numSensor ) ) ) ^ ( alpha1(i)*scaleArray(i)*numSensor - 1); 
    P1_3(i) = -P1_2(i) + alpha1(i); 
    
    Y(i) = P1_1(i) * (scaleArray(i)*numSensor + 1)/2  ...
         + P1_2(i) * (1 + 2*scaleArray(i)*numSensor + alpha1(i)*scaleArray(i)*numSensor)/2 ...
         + P1_3(i) * (scaleArray(i)*numSensor + alpha1(i)*scaleArray(i)*numSensor);
end

Y15 = zeros(1, 9);
P15_1 = zeros(1, 9);
P15_2 = zeros(1, 9);
P15_3 = zeros(1, 9);
for i = 1:9
    P15_1(i) = 1 - alpha15(i);
    P15_2(i) = alpha15(i) * (1 - ( 1 / ( alpha15(i)*scaleArray(i)*numSensor ) ) ) ^ ( alpha15(i)*scaleArray(i)*numSensor - 1); 
    P15_3(i) = -P15_2(i) + alpha15(i); 
    
    Y15(i) = P15_1(i) * (1.5*scaleArray(i)*numSensor + 1)/2  ...
         + P15_2(i) * (1 + 2*1.5*scaleArray(i)*numSensor + alpha15(i)*scaleArray(i)*numSensor)/2 ...
         + P15_3(i) * (1.5*scaleArray(i)*numSensor + alpha15(i)*scaleArray(i)*numSensor);
    %Y15(i) = (1-alpha(2))*(1.5*scaleArray(i)*numSensor + 1)/2  +  alpha(2)*(1 + 2*1.5*scaleArray(i)*numSensor + alpha(1)*scaleArray(i)*numSensor)/2;
end

Y2 = zeros(1, 9);
P2_1 = zeros(1, 9);
P2_2 = zeros(1, 9);
P2_3 = zeros(1, 9);
for i = 1:9
    P2_1(i) = 1 - alpha2(i);
    P2_2(i) = alpha2(i) * (1 - ( 1 / ( alpha2(i)*scaleArray(i)*numSensor ) ) ) ^ ( alpha2(i)*scaleArray(i)*numSensor - 1); 
    P2_3(i) = -P2_2(i) + alpha2(i); 
    
    Y2(i) = P2_1(i) * (2*scaleArray(i)*numSensor + 1)/2  ...
         + P2_2(i) * (1 + 2*2*scaleArray(i)*numSensor + alpha2(i)*scaleArray(i)*numSensor)/2 ...
         + P2_3(i) * (2*scaleArray(i)*numSensor + alpha2(i)*scaleArray(i)*numSensor);
    %Y2(i) = (1-alpha(3))*(2*scaleArray(i)*numSensor + 1)/2  +  alpha(3)*(1 + 2*2*scaleArray(i)*numSensor + alpha(3)*scaleArray(i)*numSensor)/2;
end

Y = Y*(4*10^(-3));
Y = Y+  44 * (decodingTime ) ;
Y15 = Y15*(4*10^(-3));
Y15 = Y15+  44 * (decodingTime ) ;
Y2 = Y2*(4*10^(-3));
Y2 = Y2+  44 * (decodingTime ) ;


%{
plot( scaleArray*numSensor , Transmission_Time(1,:), 'b*');
plot( scaleArray*numSensor ,Y , 'b-.');
plot( scaleArray*numSensor , Transmission_TimeR(1,:), 'b-*');

plot( scaleArray*numSensor , Transmission_Time(2,:), 'gs');
plot( scaleArray*numSensor , Y15, 'g--');
plot( scaleArray*numSensor , Transmission_TimeR(2,:), 'g-s');

plot( scaleArray*numSensor , Transmission_Time(3,:), 'rd');
plot( scaleArray*numSensor , Y2, 'r:');
plot( scaleArray*numSensor , Transmission_TimeR(3,:), 'r-d');
%}

plot( scaleArray*numSensor , Transmission_Time(1,:), 'bh');
plot( scaleArray*numSensor ,Y , 'b:', 'LineWidth', 1.2);
plot( scaleArray*numSensor , Transmission_TimeR(1,:), 'b-o');

plot( scaleArray*numSensor , Transmission_Time(2,:), 'gs');
plot( scaleArray*numSensor , Y15, 'g--');
plot( scaleArray*numSensor , Transmission_TimeR(2,:), 'g-d');

plot( scaleArray*numSensor , Transmission_Time(3,:), 'r^');
plot( scaleArray*numSensor , Y2, 'r-.');
plot( scaleArray*numSensor , Transmission_TimeR(3,:), 'r-v');




% ---------------------------------------------------------------------Reference--------------------------------------


%Y = (1- 1./scaleArray/numSensor).^(numSensor-1);
%plot( scaleArray , Y, 'b--');
%plot( scaleArray , Success, 'g:^');

lgd = legend('Ours L = N', 'Ours L = N (ana)',  'CSMA L = N', 'Ours L = 1.5N', 'Ours L = 1.5N (ana)', 'CSMA L = 1.5N', 'Ours L = 2N', 'Ours L = 2N (ana)', 'CSMA L = 2N');
lgd.Position(1) = 0.17;
lgd.Position(2) = 0.55;
lgd.FontSize = 9;

xticks(100:10:200);

set(gca,'FontSize',12);

%lgd = legend('Ours L = N','Ours L = 1.5N','Ours L =2N', 'CSMA L = N','CSMA L = 1.5N','CSMA L =2N');
%set(lgd,'Interpreter','latex');
%lgd.FontSize=12;
%lgd.FontWeight='bold';


xlabel('Number of devices');  
ylabel('Access delay (s)');

hold off;

%----------------------------------------- sub plot -------------------------------
axes('Position',[0.68 0.32 0.18 0.18]);

hold on;

plot( scaleArray*numSensor , Transmission_Time(1,:), 'bh');
plot( scaleArray*numSensor ,Y , 'b:', 'LineWidth', 1.2);

plot( scaleArray*numSensor , Transmission_Time(2,:), 'gs');
plot( scaleArray*numSensor , Y15, 'g--');

plot( scaleArray*numSensor , Transmission_Time(3,:), 'r^');
plot( scaleArray*numSensor , Y2, 'r-.');



annotation('textarrow',[0.60 0.66], [0.22 0.30]);
%annotation('arrow', [0 0], [4 4]);
xlim([160 200]);
%ylim([6*10^(-5) 8.5*10^(-5)]);
hold off;

%ylim([0 8]);

%figure1 = figure(1);
%set(figure1, 'PaperUnits', 'inches');
%set(figure1, 'PaperPosition', [0 0 4.7 5]);

print(f1, '-depsc', ['picture\','S1_TransmissionTime_Sim_', num2str(Sim_time),'_N_', num2str(num_devices),  '.eps']);
print(f1, '-dpng', ['picture\','S1_TransmissionTime_Sim_', num2str(Sim_time),'_N_', num2str(num_devices),  '.png']);