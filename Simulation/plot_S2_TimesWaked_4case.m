% Plot 
close all;
clear;
clc;
Sim_time_ppp = 400;
num_devices_ppp = 100;
load(['data/S2_simtime_', num2str(Sim_time_ppp), '_numD_', num2str(num_devices_ppp)]);

f1 = figure(1);

hold on



plot( scaleArray*num_devices , Times_Waked_prim_mod(1,:), 'r-*');
%plot( scaleArray*num_devices , Times_Waked_prim(1,:), 'r:*');
plot( scaleArray*num_devices , Times_Waked_prim_final(1,:), 'r-.*');
plot( scaleArray*num_devices , Times_Waked_convex(1,:), 'r--*');

plot( scaleArray*num_devices , Times_Waked_prim_mod(2,:), 'g-s');
%plot( scaleArray*num_devices , Times_Waked_prim(2,:), 'g:s');
plot( scaleArray*num_devices , Times_Waked_prim_final(2,:), 'g-.s');
plot( scaleArray*num_devices , Times_Waked_convex(2,:), 'g--s');

plot( scaleArray*num_devices , Times_Waked_prim_mod(3,:), 'b-d');
%plot( scaleArray*num_devices , Times_Waked_prim(3,:), 'b:d');
plot( scaleArray*num_devices , Times_Waked_prim_final(3,:), 'b-.d');
plot( scaleArray*num_devices , Times_Waked_convex(3,:), 'b--d');

%ylim([0 1]);



%Y = (1- 1./scaleArray/num_devices).^(num_devices-1);
%plot( scaleArray , Y, 'b--');
%plot( scaleArray , Success, 'g:^');
%lgd = legend('Ours L = N', 'Prim L = N', 'PrimM L = N' ,'Square L = N', 'Ours L = 1.5N', 'Prim L = 1.5N', 'PrimM L = 1.5N' , 'Square L = 1.5N', 'Ours L = 2N', 'Prim L = 2N', 'PrimM L = 2N' , 'Square L =2N');
lgd = legend('Prim L = N', 'PrimM L = N' ,'Convex L = N', 'Prim L = 1.5N', 'PrimM L = 1.5N' , 'Convex L = 1.5N', 'Prim L = 2N', 'PrimM L = 2N' , 'Convex L =2N');

% lgd = legend('Ours L = N','Ours L = 1.5N','Ours L = 2N', 'Ours (ana) L = N \alpha = 0.49', 'Ours (ans) L = 1.5N \alpha = 0.34', 'Ours (ana) L = 2N, \alpha = 0.24', 'CSMA L = N','CSMA L = 1.5N','CSMA L =2N');
%set(lgd,'Interpreter','latex');
%lgd.FontSize=12;
%lgd.FontWeight='bold';


xlabel('Number of devices');  
ylabel('Times of being waked');

hold off;

print(f1, '-depsc', ['picture\','S2_4case_TimesWaked_Sim_', num2str(Sim_time_ppp),'_N_', num2str(num_devices_ppp),  '.eps']);