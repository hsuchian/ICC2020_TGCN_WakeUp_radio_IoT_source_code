% Plot 
close all;
clear;
clc;
Sim_time_ppp = 1800;
num_devices_ppp = 100;
load(['data/S2_simtime_', num2str(Sim_time_ppp), '_numD_', num2str(num_devices_ppp)]);

f1 = figure(1);
hold on;


plot( scaleArray*num_devices , Transmission_Time_prim_final(1,:), 'b-o');
plot( scaleArray*num_devices , Transmission_Time_prim_mod(1,:), 'b--o');
plot( scaleArray*num_devices , Transmission_Time_sq(1,:), 'b-.o');

plot( scaleArray*num_devices , Transmission_Time_prim_final(2,:), 'g-s');
plot( scaleArray*num_devices , Transmission_Time_prim_mod(2,:), 'g--s');
plot( scaleArray*num_devices , Transmission_Time_sq(2,:), 'g-.s');

plot( scaleArray*num_devices , Transmission_Time_prim_final(3,:), 'r-d');
plot( scaleArray*num_devices , Transmission_Time_prim_mod(3,:), 'r--d');
plot( scaleArray*num_devices , Transmission_Time_sq(3,:), 'r-.d');

%ylim([0 5]);



%Y = (1- 1./scaleArray/num_devices).^(num_devices-1);
%plot( scaleArray , Y, 'b--');
%plot( scaleArray , Success, 'g:^');

lgd = legend('ALG2 L = N', 'ALG1 L = N', 'Square L = N', 'ALG2 L = 1.5N', 'ALG1 L = 1.5N', 'Square L = 1.5N', 'ALG2 L = 2N', 'ALG1 L = 2N', 'Square L = 2N');
lgd.Position(1) = 0.17;
lgd.Position(2) = 0.53;

lgd.FontSize=9;
set(gca,'FontSize',12);

xticks(200:50:700);

%lgd.FontWeight='bold';


xlabel('Number of devices');  
ylabel('Access delay (s)');

hold off;


print(f1, '-depsc', ['picture\','S2_TransmissionTime_Sim_', num2str(Sim_time_ppp),'_N_', num2str(num_devices_ppp),  '.eps']);
print(f1, '-dpng', ['picture\','S2_TransmissionTime_Sim_', num2str(Sim_time),'_N_', num2str(num_devices),  '.png']);