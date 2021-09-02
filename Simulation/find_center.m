
function [mid_x, mid_y] = find_center(X,Y, UAVradius)

    x0 = [mean(X), mean(Y)];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    
    objectivef = @(x) 1;
    nlcon1 = @(x) nclon(x, X,Y,UAVradius);
    options = optimoptions('fmincon', 'Display','off', 'MaxIterations', 10000);
    
    [x_sol, fval, exitflag] = fmincon(objectivef,x0,A,b,Aeq,beq,lb,ub, nlcon1, options); 
    
    if exitflag == 0
        fprintf('Exceed Tolerance\n'); 
    end
    
    mid_x = x_sol(1);
    mid_y = x_sol(2);
    
    %{
    close all;
    figure(1)
    hold on;
    plot(X, Y, 'r.');
    
    for i = 1:length(X)
        circle(UAVradius, X(i), Y(i));
    end
    plot(mid_x, mid_y, 'b.');
    hold off
    %}
end