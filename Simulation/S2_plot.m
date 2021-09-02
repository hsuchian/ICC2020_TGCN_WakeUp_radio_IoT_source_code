% plot even partition

hold on
plot(pointX, pointY, '.');

pointX = x;
pointY= y;

width = 20*2^(1/2);
for i = 1:100/width
    plot(pointX, width*i*ones(1,length(pointY)), 'r-');
end


for i = 1:100/width
    plot(width*i*ones(1,length(pointX)),pointY, 'r-');
end

plot(midpoint(1,:), midpoint(2,:), 'g.')