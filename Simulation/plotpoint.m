
hold on
x = pointX;
y=pointY;
plot(x,y, '.')

 for i = 1:20
     text(x(i)+0.5, y(i)+0.5, num2str(i))
 end