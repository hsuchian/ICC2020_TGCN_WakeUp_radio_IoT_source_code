

% nlcon1 = @(x) nlcon(x, m, n, P_constraint, N, Current_load, length_subframe);

function [c, ceq] = nclon(x, X, Y, UAVradius)
    ceq = [];
    
    num_X = length(X);
    c = zeros(1, num_X);
    
    for i = 1:num_X
        c(i) = (X(i) - x(1))^2  +  (Y(i) - x(2))^2  -  UAVradius^2;
    end
end