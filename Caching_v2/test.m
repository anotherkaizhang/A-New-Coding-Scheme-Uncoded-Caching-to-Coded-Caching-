Ned = min(K,N); % only consider demand types where all files are requested
figure
% Tian
t_axis = 0:K;
M_Tian = t_axis.*((N-1)*t_axis+K-N)/K/(K-1);
R_Tian =  N*(K-t_axis)/K;
plot(M_Tian,R_Tian, 'm-*', 'Linewidth', 2) 

% Yu 
hold on
t_axis = 0:K-1;
M_Yu = t_axis*N/K;
R_Yu = zeros(size(t_axis)); 
for iR = 1:length(t_axis)
    if K - Ned >= t_axis(iR)+1
        R_Yu(iR) = (nchoosek(K, t_axis(iR)+1) - nchoosek(K - Ned, t_axis(iR)+1)) / nchoosek(K,t_axis(iR));
    else
        R_Yu(iR) = nchoosek(K, t_axis(iR)+1) / nchoosek(K, t_axis(iR));
    end
end
M_Yu = [M_Yu, N];
R_Yu = [R_Yu, 0];
plot(M_Yu,R_Yu,'b-o', 'Linewidth', 2) 
legend('Tian','Yu')
% =====================================

%% t = 2 black line:
% M = [15, 11]/15;
% R = [19, 25]/15;
% plot(M,R,'g-d', 'Linewidth', 1.1) 
% 
% M = [15, 11]/15;
% R = [19, 24]/15;
% plot(M,R,'g-d', 'Linewidth', 1.1) 

%% t = 3 green line:
% M = [30, 28]/20;
% R = [15, 18]/20;
% plot(M,R,'g-d', 'Linewidth', 1.1) 
% 
% M = [30, 20]/20;
% R = [15, 27]/20;
% plot(M,R,'g-d', 'Linewidth', 1.1) 
% =====================================
load('P_all.mat')
temp = P_all{3,1}; temp = [temp; temp(1,:)]; x1 = temp(:,1), y1 = temp(:,2);
temp = P_all{3,2}; temp = [temp; temp(1,:)]; x2 = temp(:,1), y2 = temp(:,2);
temp = P_all{3,3}; temp = [temp; temp(1,:)]; x3 = temp(:,1), y3 = temp(:,2);

[xb, yb] = polybool('intersection', x1, y1, x2, y2);
[xb, yb] = polybool('intersection', xb, yb, x3, y3);
xb(end) = []; yb(end) = [];
plot(xb, yb, 'g-d')
% % patch(xb, yb, 1, 'FaceColor', 'r')
% axis equal, axis off, hold on
% % plot(x1, y1, x2, y2, x3, y3, 'Color', 'k')
% title('Intersection')


% [1,1.26666666666667;0.933333333333333,1.33333333333333;0.466666666666667,2]
% x1 = cos(theta) - 0.5;
% y1 = -sin(theta);    % -sin(theta) to make a clockwise contour
% x2 = x1 + 1;
% y2 = y1;
% [xa, ya] = polybool('union', x1, y1, x2, y2);
% [xb, yb] = polybool('intersection', x1, y1, x2, y2);
% [xc, yc] = polybool('xor', x1, y1, x2, y2);
% [xd, yd] = polybool('subtraction', x1, y1, x2, y2);
% 
% subplot(2, 2, 1)
% patch(xa, ya, 1, 'FaceColor', 'r')
% axis equal, axis off, hold on
% plot(x1, y1, x2, y2, 'Color', 'k')
% title('Union')
% 
% subplot(2, 2, 2)
% patch(xb, yb, 1, 'FaceColor', 'r')
% axis equal, axis off, hold on
% plot(x1, y1, x2, y2, 'Color', 'k')
% title('Intersection')
% 
% subplot(2, 2, 3)
% % The output of the exclusive-or operation consists of disjoint
% % regions.  It can be plotted as a single patch object using the
% % face-vertex form.  Use poly2fv to convert a polygonal region
% % to face-vertex form.
% [f, v] = poly2fv(xc, yc);
% patch('Faces', f, 'Vertices', v, 'FaceColor', 'r', ...
%   'EdgeColor', 'none')
% axis equal, axis off, hold on
% plot(x1, y1, x2, y2, 'Color', 'k')
% title('Exclusive Or')
% 
% subplot(2, 2, 4)
% patch(xd, yd, 1, 'FaceColor', 'r')
% axis equal, axis off, hold on
% plot(x1, y1, x2, y2, 'Color', 'k')
% title('Subtraction')