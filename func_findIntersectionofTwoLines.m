function point_intersect = func_findIntersectionofTwoLines(point1, point2, point3, point4)
%line1
% x1  = [7.8 8.5];
x1 = [point1(1), point2(1)];
% y1  = [0.96 0.94];
y1 = [point1(2), point2(2)];
%line2
x2 = [point3(1), point4(1)];
% x2 = [8.25 8.25];
y2 = [point3(2), point4(2)];
% y2 = [0 0.99];
%fit linear polynomial
p1 = polyfit(x1,y1,1);
p2 = polyfit(x2,y2,1);
%calculate intersection
x_intersect = fzero(@(x) polyval(p1-p2,x),3);
y_intersect = polyval(p1,x_intersect);
point_intersect = [x_intersect, y_intersect];
% figure
% line(x1,y1);
% hold on;
% line(x2,y2);
% plot(x_intersect,y_intersect,'r*')