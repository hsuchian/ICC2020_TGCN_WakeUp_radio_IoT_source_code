% plot


PLSleep  = 1.9*(10^(-6));
PActive = 4.1*(10^(-3));

k = 0.8;       % percentage hashing collision
N = 100;
m = 100;
L = 50:5:200;     %
Time = zeros(1, length(L));
Energy = zeros(1, length(L));
SuccessProb = zeros(1, length(L));
for i = 1:length(L)
    Time(i) = (1-k)*(1+L(i))/2 + k*(1+2*L(i)+m)/2;
    Energy(i) = (1-k)*(PActive + (L(i)-1)*PLSleep) + k*(2*PActive + PLSleep * (  L(i) - 2 + (1+m)/2)  );
    SuccessProb(i) = (1-k) + k*(  (1- (1/m))^(N*k-1) );
end



plot(L, Time)
figure(2)
plot(L, Energy, 'r')
figure(3)
plot(L, SuccessProb, 'g')