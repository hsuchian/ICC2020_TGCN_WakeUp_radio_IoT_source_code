
clear ;
clc;

temp= 6;
num_device =temp;
num_slot = temp;
syms x
prob_n = zeros(1, num_device);
prob_cum_n = zeros(1, num_device);

for k = 1:3
    if (k ~= num_device)
        prob_n(k) =  nchoosek(num_device, k) * nchoosek(num_slot, k) * factorial(k) * (num_slot - k)^(num_device-k);
        fprintf('probn_1 = %d\n', prob_n(k));
        prob_n_next(k) =  (nchoosek(num_device, k+1) * nchoosek(num_slot, k+1) * factorial(k+1) * (num_slot-k-1)^(num_device-k-1));
        fprintf('probn_2 = %d\n', prob_n_next(k));
        %prob_n(k) = prob_n(k) / ((num_slot)^(num_device));
    else 
        prob_n(k) = factorial(k);
        %prob_n(k) = prob_n(k) / ((num_slot)^(num_device));
    end
    prob_cum_n(k) = sum(prob_n(1:k));
end


plot(1:temp, prob_n(1:temp));
%plot(1:10, prob_cum_n(1:10));
