% x = 30:10:100;
% figure();
% % axes('NextPlot', 'add');  % Same as: hold on
% hold on
% for k = 1:8
%     y{k} = rand(size(x));
%     plot(x, y{k}, 'Color', rand(1, 3))
%     %axis equal
% end
x1=[2 3 4 5];
y1=[9 4 3 2];
x2=[11 20 30 50 ];
y2= [ 20 30 50 60];
figure()
axis equal
plot(x1,y1) 
hold on 
plot(x2,y2)
