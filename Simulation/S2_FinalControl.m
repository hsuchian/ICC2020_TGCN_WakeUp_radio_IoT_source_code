

% S2_Control
clear;
clc;

global hash1
hash1 = [1 1 1 1 1 1 1            8.636364e-01 7.979798e-01 8.000000e-01 7.851240e-01 7.878788e-01 7.622378e-01 7.272727e-01 6.969697e-01 7.215909e-01 6.577540e-01 6.969697e-01 6.889952e-01 6.636364e-01 6.709957e-01 6.487603e-01 6.284585e-01 6.325758e-01 6.581818e-01 6.153846e-01 6.161616e-01 6.071429e-01 6.206897e-01 6.121212e-01 6.011730e-01 6.022727e-01 5.840220e-01 6.016043e-01 5.818182e-01 5.808081e-01 5.724816e-01 5.598086e-01 5.780886e-01 5.772727e-01 5.631929e-01 5.822511e-01 5.708245e-01 5.888430e-01 5.555556e-01 5.612648e-01 5.415861e-01 5.625000e-01 5.602968e-01 5.509091e-01 5.561497e-01 5.629371e-01 5.506003e-01 5.555556e-01 5.553719e-01 5.503247e-01 5.454545e-01 5.517241e-01 5.315871e-01 5.303030e-01 5.380030e-01 5.263930e-01 5.209235e-01 5.241477e-01 5.384615e-01 5.344353e-01 5.345997e-01 5.320856e-01 5.401845e-01 5.220779e-01 5.185659e-01 5.214646e-01 5.242839e-01 5.257985e-01 5.200000e-01 5.215311e-01 5.064935e-01 5.069930e-01 5.258918e-01 5.193182e-01 5.072952e-01 5.133038e-01 5.027382e-01 5.054113e-01 5.315508e-01 5.084567e-01 5.078370e-01 5.030992e-01 5.035751e-01 5.000000e-01 5.014985e-01 5.128458e-01 5.083089e-01 5.038685e-01 5.023923e-01 5.075758e-01 4.985942e-01 5.064935e-01 5.050505e-01 4.918182e-01 4.932493e-01 4.964349e-01 4.986761e-01 5.017483e-01 4.935065e-01 4.939966e-01 4.842821e-01 4.957912e-01 5.012510e-01 4.900826e-01 4.889435e-01 4.918831e-01 4.883347e-01 4.928230e-01 4.837945e-01 4.874608e-01 4.879565e-01 4.861325e-01 4.843392e-01 5.007576e-01        4.848485e-01 4.881603e-01 4.832882e-01 4.847670e-01 4.782222e-01 4.753086e-01 4.899388e-01 4.739583e-01 4.728682e-01 4.752137e-01 4.826124e-01 4.949495e-01 4.703425e-01 4.809287e-01 4.962963e-01 4.803922e-01 4.768856e-01 4.718196e-01 4.796163e-01 4.809524e-01 4.783294e-01 4.757433e-01 4.809635e-01 4.791667e-01 4.781609e-01 4.779300e-01 4.746788e-01 4.804805e-01 4.683072e-01 4.666667e-01 4.672553e-01 4.736842e-01 4.684096e-01 4.668110e-01 4.788530e-01 4.750712e-01 4.713376e-01 4.873418e-01 4.744934e-01 4.729167e-01 4.720497e-01 4.643347e-01 4.635310e-01 4.701897e-01 4.929293e-01 4.558233e-01 4.763806e-01 4.576720e-01 4.687705e-01 4.679739e-01 4.671865e-01 4.625323e-01 4.669236e-01 4.610473e-01 4.653968e-01 4.671717e-01 4.607659e-01 4.719101e-01 4.655493e-01 4.641975e-01 4.647023e-01 4.743590e-01 4.632665e-01 4.553140e-01 4.660661e-01 4.629630e-01 4.670232e-01 4.592199e-01 4.697237e-01 4.637427e-01 4.624782e-01 4.646991e-01 4.565343e-01 4.599084e-01 4.592593e-01 4.614512e-01 4.692611e-01 4.612795e-01 4.595198e-01 4.661111e-01 ]; 
global hash15
hash15 = [1 1 1 1 1 1 1 1 1 1 1 9.090909e-01 8.741259e-01 9.350649e-01 8.787879e-01 9.090909e-01 8.930481e-01 8.989899e-01 8.516746e-01 8.454545e-01 8.225108e-01 8.429752e-01 8.102767e-01 8.068182e-01 8.145455e-01 7.972028e-01 8.148148e-01 8.181818e-01 7.742947e-01 7.696970e-01 7.859238e-01 7.642045e-01 7.768595e-01 7.673797e-01 7.792208e-01 7.676768e-01 7.690418e-01 7.559809e-01 7.482517e-01 7.454545e-01 7.427938e-01 7.402597e-01 7.547569e-01 7.500000e-01 7.313131e-01 7.391304e-01 7.427466e-01 7.518939e-01 7.291280e-01 7.218182e-01 7.326203e-01 7.202797e-01 7.409949e-01 7.222222e-01 6.975207e-01 7.029221e-01 7.081340e-01 6.927900e-01 7.241911e-01 6.909091e-01 7.093890e-01 7.082111e-01 6.868687e-01 7.215909e-01 6.867133e-01 6.997245e-01 7.028494e-01 6.778075e-01 6.877470e-01 6.909091e-01 6.901408e-01 6.969697e-01 6.861768e-01 6.842752e-01 6.884848e-01 6.746411e-01 6.871311e-01 6.876457e-01 6.823936e-01 6.784091e-01 6.891134e-01 6.629712e-01 6.648412e-01 6.807359e-01 6.695187e-01 6.606765e-01 6.781609e-01 6.725207e-01 6.772217e-01 6.939394e-01 6.703297e-01 6.689723e-01 6.627566e-01 6.634429e-01 6.574163e-01 6.742424e-01 6.569822e-01 6.697588e-01 6.565657e-01 6.636364e-01 6.579658e-01 6.604278e-01 6.548985e-01 6.625874e-01 6.597403e-01 6.586621e-01 6.499575e-01 6.675084e-01 6.480400e-01 6.595041e-01 6.486486e-01 6.534091e-01 6.452132e-01 6.435407e-01 6.505929e-01 6.528213e-01 6.418026e-01 6.563945e-01 6.524064e-01 6.446970e-01                       6.473829e-01 6.502732e-01 6.603433e-01 6.550179e-01 6.346667e-01 6.490300e-01 6.430446e-01 6.614583e-01 6.339363e-01 6.478632e-01 6.649703e-01 6.414141e-01 6.574770e-01 6.451078e-01 6.436214e-01 6.429739e-01 6.496350e-01 6.417069e-01 6.442846e-01 6.460317e-01 6.390859e-01 6.377152e-01 6.472416e-01 6.373457e-01 6.398467e-01 6.400304e-01 6.371882e-01 6.381381e-01 6.390753e-01 6.407407e-01 6.320824e-01 6.352339e-01 6.499637e-01 6.327561e-01 6.243728e-01 6.360399e-01 6.433121e-01 6.230661e-01 6.429071e-01 6.250000e-01 6.287095e-01 6.474623e-01 6.264485e-01 6.300813e-01 6.329966e-01 6.258367e-01 6.340652e-01 6.203704e-01 6.344510e-01 6.294118e-01 6.452242e-01 6.259690e-01 6.152858e-01 6.321839e-01 6.330159e-01 6.363636e-01 6.283741e-01 6.317104e-01 6.244569e-01 6.172840e-01 6.310620e-01 6.172161e-01 6.290225e-01 6.195652e-01 6.264264e-01 6.302270e-01 6.268568e-01 6.223404e-01 6.243386e-01 6.233918e-01 6.352531e-01 6.261574e-01 6.275187e-01 6.231386e-01 6.216524e-01 6.184807e-01 6.277496e-01 6.167228e-01 6.203238e-01 6.255556e-01 ];
global hash2
hash2 =  [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1            1 9.251337e-01 9.696970e-01 9.330144e-01 9.090909e-01 9.134199e-01 9.173554e-01 9.209486e-01 9.090909e-01 9.054545e-01 9.230769e-01 8.653199e-01 8.896104e-01 8.714734e-01 8.848485e-01 9.002933e-01 8.750000e-01 8.429752e-01 8.449198e-01 8.545455e-01 8.737374e-01 8.648649e-01 8.540670e-01 8.414918e-01 8.727273e-01 8.447894e-01 8.419913e-01 8.350951e-01 8.326446e-01 8.303030e-01 8.122530e-01 8.355899e-01 8.125000e-01 8.126160e-01 8.236364e-01 8.074866e-01 8.234266e-01 8.096055e-01 8.148148e-01 8.049587e-01 8.214286e-01 8.006380e-01 8.056426e-01 7.996918e-01 8.121212e-01 8.062593e-01 7.815249e-01 7.979798e-01 8.039773e-01 7.902098e-01 7.920110e-01 7.978290e-01 7.941176e-01 7.667984e-01 7.935065e-01 7.861716e-01 7.878788e-01 7.907846e-01 7.825553e-01 7.733333e-01 7.906699e-01 7.744982e-01 7.855478e-01 7.733026e-01 7.795455e-01 7.811448e-01 7.716186e-01 7.667032e-01 7.824675e-01 7.732620e-01 7.758985e-01 7.763845e-01 7.685950e-01 7.732380e-01 7.656566e-01 7.672328e-01 7.638340e-01 7.644184e-01 7.427466e-01 7.540670e-01 7.623106e-01 7.666354e-01 7.560297e-01 7.548209e-01 7.609091e-01 7.596760e-01 7.575758e-01 7.502207e-01 7.613636e-01 7.549784e-01 7.444254e-01 7.561597e-01 7.550505e-01 7.597998e-01 7.520661e-01 7.452907e-01 7.556818e-01 7.530169e-01 7.408293e-01 7.494071e-01 7.586207e-01 7.482517e-01 7.488444e-01 7.540107e-01 7.272727e-01             7.483930e-01 7.440801e-01 7.262873e-01 7.419355e-01 7.413333e-01 7.425044e-01 7.489064e-01 7.413194e-01 7.338501e-01 7.316239e-01 7.599661e-01 7.449495e-01 7.468672e-01 7.330017e-01 7.374486e-01 7.320261e-01 7.485807e-01 7.286634e-01 7.306155e-01 7.452381e-01 7.375887e-01 7.425665e-01 7.428127e-01 7.368827e-01 7.409962e-01 7.420091e-01 7.294029e-01 7.409910e-01 7.330350e-01 7.488889e-01 7.343635e-01 7.295322e-01 7.349310e-01 7.236652e-01 7.369176e-01 7.364672e-01 7.211607e-01 7.447257e-01 7.267645e-01 7.222222e-01 7.149758e-01 7.229081e-01 7.191547e-01 7.256098e-01 7.272727e-01 7.275770e-01 7.238856e-01 7.321429e-01 7.186062e-01 7.248366e-01 7.244964e-01 7.273902e-01 7.174053e-01 7.113665e-01 7.161905e-01 7.070707e-01 7.150031e-01 7.253433e-01 7.219119e-01 7.277778e-01 7.237569e-01 7.130647e-01 7.255616e-01 7.131643e-01 7.171171e-01 7.138590e-01 7.254902e-01 7.151300e-01 7.189888e-01 7.257310e-01 7.219314e-01 7.193287e-01 7.213587e-01 7.222222e-01 7.219373e-01 7.233560e-01 7.241963e-01 7.098765e-01 7.236181e-01 7.161111e-01 ];


Sim_time = 1800;

num_devices = 100;
max_range = 200;       % length width in meter
UAVradius = 20;        % meter
width = UAVradius*(2^(0.5));

% parameters
scaleArray = 2:0.25:7;
scaleArray_frame_length = 1:0.5:2;
Length_of_hash = length(hash2);

tic 
output_length = (7-2)/0.25 +1;

Success_Prob_prim_mod = zeros(3, output_length ); 
Transmission_Time_prim_mod = zeros(3, output_length ); 
Energy_Consumption_prim_mod = zeros(3, output_length ); 
Times_Waked_prim_mod = zeros(3, output_length ); 


Success_Prob_prim_final = zeros(3, output_length ); 
Transmission_Time_prim_final = zeros(3, output_length ); 
Energy_Consumption_prim_final = zeros(3, output_length ); 
Times_Waked_prim_final = zeros(3, output_length ); 

Success_Prob_convex = zeros(3, output_length ); 
Transmission_Time_convex = zeros(3, output_length ); 
Energy_Consumption_convex = zeros(3, output_length ); 
Times_Waked_convex = zeros(3, output_length ); 

Success_Prob_sq = zeros(3, output_length ); 
Transmission_Time_sq = zeros(3, output_length ); 
Energy_Consumption_sq = zeros(3, output_length ); 
Times_Waked_sq = zeros(3, output_length ); 

% start simulation
index1 = 1;       % frame length sclae
index2 = 1;       % devices number scale

 for scale = scaleArray 
     for i = 1:Sim_time
        
        % generate random points
        pointX = rand(1, round(num_devices*scale)*10);
        pointY = rand(1, round(num_devices*scale)*10);
        pointX = pointX*max_range;
        pointY = pointY*max_range;
            
            
        % generate shape         
        x_gen = rand(1,30)*max_range;
        y_gen = rand(1,30)*max_range;
        
        shape = polyshape(x_gen,y_gen);
        
        in_or_not = isinterior(shape, pointX, pointY);
        in_or_not_int = transpose(int8(in_or_not));

        feasible_ind = find(in_or_not_int == 1);
        
        pointX = pointX(feasible_ind);   
        pointY = pointY(feasible_ind); 
        if length(feasible_ind) > round(num_devices*scale)
            pointX = pointX(1:round(num_devices*scale));
            pointY = pointY(1:round(num_devices*scale));
        else 
            fprintf('point not enough %d %d\n',length(feasible_ind) ,round(num_devices*scale) );
            pause
        end
        % check if all random points reside in polygon
        %for id = 1:length(pointY)
        %    while ~(isinterior(shape,pointX(id),pointY(id)))
        %        pointX(id) = rand()*max_range;
        %        pointY(id) = rand()*max_range;
        %    end
        %end
            
            
        index1 = 1;   % frame length sclae index   
        for fscale = scaleArray_frame_length 
            % S2_mainfunction(max_range, num_devices, UAVradius, fscale, pointX, pointY)
        
            [SuccessProb_prim_mod, TransmissionTime_prim_mod, EnergyConsumption_prim_mod, TimesWaked_prim_mod ] = ...
            S2_mainfunction_prim_mod(max_range, round(num_devices*scale), UAVradius, fscale, pointX, pointY, Length_of_hash);
        
            [SuccessProb_prim_final, TransmissionTime_prim_final, EnergyConsumption_prim_final, TimesWaked_prim_final ] = ...
            S2_mainfunction_prim_final(max_range, round(num_devices*scale), UAVradius, fscale, pointX, pointY, Length_of_hash);
        
            %[SuccessProb_convex, TransmissionTime_convex, EnergyConsumption_convex, TimesWaked_convex ] = ...
            %S2_mainfunction_convex(max_range, round(num_devices*scale), UAVradius, fscale, pointX, pointY, Length_of_hash);
        
            [SuccessProb_sq, TransmissionTime_sq, EnergyConsumption_sq, TimesWaked_sq] = ...
            S2_mainfunction_sq(max_range, round(num_devices*scale), UAVradius, fscale, pointX, pointY, Length_of_hash);
            
            
            Success_Prob_prim_mod(index1, index2)       = Success_Prob_prim_mod(index1, index2) + SuccessProb_prim_mod;
            Transmission_Time_prim_mod(index1, index2)  = Transmission_Time_prim_mod(index1, index2) + TransmissionTime_prim_mod;
            Energy_Consumption_prim_mod(index1, index2) = Energy_Consumption_prim_mod(index1, index2) + EnergyConsumption_prim_mod;
            Times_Waked_prim_mod(index1, index2)        = Times_Waked_prim_mod(index1, index2) + TimesWaked_prim_mod; 
            
            Success_Prob_prim_final(index1, index2)       = Success_Prob_prim_final(index1, index2) + SuccessProb_prim_final;
            Transmission_Time_prim_final(index1, index2)  = Transmission_Time_prim_final(index1, index2) + TransmissionTime_prim_final;
            Energy_Consumption_prim_final(index1, index2) = Energy_Consumption_prim_final(index1, index2) + EnergyConsumption_prim_final;
            Times_Waked_prim_final(index1, index2)        = Times_Waked_prim_final(index1, index2) + TimesWaked_prim_final; 
            
            %Success_Prob_convex(index1, index2)       = Success_Prob_convex(index1, index2) + SuccessProb_convex;
            %Transmission_Time_convex(index1, index2)  = Transmission_Time_convex(index1, index2) + TransmissionTime_convex;
            %Energy_Consumption_convex(index1, index2) = Energy_Consumption_convex(index1, index2) + EnergyConsumption_convex;
            %Times_Waked_convex(index1, index2)        = Times_Waked_convex(index1, index2) + TimesWaked_convex; 
            
            Success_Prob_sq(index1, index2)       = Success_Prob_sq(index1, index2) + SuccessProb_sq;
            Transmission_Time_sq(index1, index2)  = Transmission_Time_sq(index1, index2) + TransmissionTime_sq;
            Energy_Consumption_sq(index1, index2) = Energy_Consumption_sq(index1, index2) + EnergyConsumption_sq;
            Times_Waked_sq(index1, index2)        = Times_Waked_sq(index1, index2) + TimesWaked_sq; 
            
            
            index1 = index1 +1; 
        end 
     end
    index2 = index2 + 1;
 end
  
  Success_Prob_prim_mod       = Success_Prob_prim_mod/Sim_time;
  Transmission_Time_prim_mod  = Transmission_Time_prim_mod/Sim_time;
  Energy_Consumption_prim_mod = Energy_Consumption_prim_mod/Sim_time;
  Times_Waked_prim_mod        = Times_Waked_prim_mod/Sim_time; 
  
  Success_Prob_prim_final       = Success_Prob_prim_final/Sim_time;
  Transmission_Time_prim_final  = Transmission_Time_prim_final/Sim_time;
  Energy_Consumption_prim_final = Energy_Consumption_prim_final/Sim_time;
  Times_Waked_prim_final        = Times_Waked_prim_final/Sim_time; 
  
  %Success_Prob_convex       = Success_Prob_convex/Sim_time;
  %Transmission_Time_convex  = Transmission_Time_convex/Sim_time;
  %Energy_Consumption_convex = Energy_Consumption_convex/Sim_time;
  %Times_Waked_convex        = Times_Waked_convex/Sim_time; 
  
  Success_Prob_sq        = Success_Prob_sq/Sim_time;
  Transmission_Time_sq   = Transmission_Time_sq/Sim_time;
  Energy_Consumption_sq  = Energy_Consumption_sq/Sim_time;
  Times_Waked_sq         = Times_Waked_sq/Sim_time;
 
toc

save(['data/S2_simtime_', num2str(Sim_time), '_numD_', num2str(num_devices)]);
%[SuccessProb, TransmissionTime, EnergyConsumption,] = S2_mainfunction(max_range, num_devices, UAVradius, fscale, pointX, pointY);
%[SuccessProb_sq, TransmissionTime_sq, EnergyConsumption_sq] = S2_mainfunction_sq(max_range, num_devices, UAVradius, fscale, pointX, pointY);
