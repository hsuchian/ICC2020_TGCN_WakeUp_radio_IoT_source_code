% Plot 
clear;
clc;
close all;

Sim_time = 150;
num_devices = 50;
load(['data/S1_simtime_', num2str(Sim_time), '_numD_', num2str(num_devices)]);

f1 = figure(1);
hold on

%hash2(150) = 0.7380;



%plot( scaleArray*numSensor , Success_Prob(2,:), 'g*');
%plot( scaleArray*numSensor , Success_Prob(3,:), 'b*');

alpha1 = 1 - hash1;
alpha1 = alpha1(round(scaleArray*numSensor));

alpha15 = 1 - hash15;
alpha15 = alpha15(round(scaleArray*numSensor));

alpha2 = 1 - hash2;
alpha2 = alpha2(round(scaleArray*numSensor));

Y = zeros(1, 9);
for i = 1:9
    Y(i) = (1-alpha1(i)) + alpha1(i)*(1 - ( 1 / (alpha1(i)*scaleArray(i)*numSensor) ))^(alpha1(i)*scaleArray(i)*numSensor -1);
end

Y15 = zeros(1, 9);
for i = 1:9
    Y15(i) = (1-alpha15(i)) + alpha15(i)*(1 - ( 1 / (alpha15(i)*scaleArray(i)*numSensor) ))^(alpha15(i)*scaleArray(i)*numSensor -1);
end

Y2 = zeros(1, 9);
for i = 1:9
    Y2(i) = (1-alpha2(i)) + alpha2(i)*(1 - ( 1 / (alpha2(i)*scaleArray(i)*numSensor) ))^(alpha2(i)*scaleArray(i)*numSensor -1);
end

%{
plot( scaleArray*numSensor , Success_Prob(1,:), 'b*');
plot( scaleArray*numSensor ,Y , 'b-.');
plot( scaleArray*numSensor , Success_ProbR(1,:), 'b-*');

plot( scaleArray*numSensor , Success_Prob(2,:), 'gs');
plot( scaleArray*numSensor , Y15, 'g--');
plot( scaleArray*numSensor , Success_ProbR(2,:), 'g-s');

plot( scaleArray*numSensor , Success_Prob(3,:), 'rd');
plot( scaleArray*numSensor , Y2, 'r:');
plot( scaleArray*numSensor , Success_ProbR(3,:), 'r-d');
%}


plot( scaleArray*numSensor , Success_Prob(1,:), 'bh');
plot( scaleArray*numSensor ,Y , 'b:', 'LineWidth', 1.2);
plot( scaleArray*numSensor , Success_ProbR(1,:), 'b-o');

plot( scaleArray*numSensor , Success_Prob(2,:), 'gs');
plot( scaleArray*numSensor , Y15, 'g--');
plot( scaleArray*numSensor , Success_ProbR(2,:), 'g-d');

plot( scaleArray*numSensor , Success_Prob(3,:), 'r^');
plot( scaleArray*numSensor , Y2, 'r-.');
plot( scaleArray*numSensor , Success_ProbR(3,:), 'r-v');


%plot( scaleArray*numSensor , Success_ProbR(1,:), 'r--*');
%plot( scaleArray*numSensor , Success_ProbR(2,:), 'g--s');


%Y = (1- 1./scaleArray/numSensor).^(numSensor-1);
%plot( scaleArray , Y, 'b--');
%plot( scaleArray , Success, 'g:^');
lgd = legend('Ours L = N', 'Ours L = N (ana)',  'CSMA L = N', 'Ours L = 1.5N', 'Ours L = 1.5N (ana)', 'CSMA L = 1.5N', 'Ours L = 2N', 'Ours L = 2N (ana)', 'CSMA L = 2N');

% lgd = legend('Ours L = N','Ours L = 1.5N','Ours L = 2N', 'Ours (ana) L = N \alpha = 0.49', 'Ours (ans) L = 1.5N \alpha = 0.34', 'Ours (ana) L = 2N, \alpha = 0.24', 'CSMA L = N','CSMA L = 1.5N','CSMA L =2N');
%set(lgd,'Interpreter','latex');
%lgd.FontSize=12;
%lgd.FontWeight='bold';

lgd.Position(1) = 0.57;
lgd.Position(2) = 0.27;
lgd.FontSize = 9;

xticks(100:10:200);
yticks(0: 0.1:0.9);

xlabel('Number of devices');  
ylabel('Success probability');

% set x-axis and y-axis font size to 14
set(gca,'FontSize',12);



hold off;
%ylim([0 1]);

%figure1 = figure(1);
%set(figure1, 'PaperUnits', 'inches');
%set(figure1, 'PaperPosition', [0 0 4.7 5]);

print(f1, '-depsc', ['picture\','S1_SuccessProb_Sim_', num2str(Sim_time),'_N_', num2str(num_devices),  '.eps']);
print(f1, '-dpng', ['picture\','S1_SuccessProb_Sim_', num2str(Sim_time),'_N_', num2str(num_devices),  '.png']);