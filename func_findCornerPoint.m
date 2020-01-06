function result = func_findCornerPoint(N, num_var_total, A,b,Aeq,beq, options, corner_points, ifAlreadyTeated)
pos = find(ifAlreadyTeated == 0, 1);
while ~isempty(pos)
    left_point = corner_points(pos, :);
    right_point = corner_points(pos+1,:);
    f = [zeros(num_var_total, 1); (right_point(2) - left_point(2))/(left_point(1) - right_point(1)); 1];
    lb = [zeros(num_var_total, 1); 0; 0]; % left_point(1); right_point(2)];
    ub = [ones(num_var_total, 1); N; N]; % right_point(1); left_point(2)];
    x = linprog(f,sparse(A),b,sparse(Aeq), beq,lb, ub, options);
%     save('x.mat', 'x');
%     figure
%     hold on
%     plot(x)
%     legend
    M = x(end-1);
    R = x(end);
%     if (abs(M - 0.9333) + abs(R-1.9467) < 0.001) % || 
%         if (abs(M-1.067)+abs(R-1.7333) < 0.001)
%         save('x.mat', 'x');
%         find(x>0)
%         x(find(x>0))
%         pause
%     end
    if func_point_to_line_distance([M,R,0], [left_point, 0], [right_point, 0]) > 1e-6 % neither left point nor right point
        corner_points = [corner_points(1:pos,:); [M,R]; corner_points(pos+1:end,:)]; 
        ifAlreadyTeated = [ifAlreadyTeated(1:pos-1), 0, 0, ifAlreadyTeated(pos+1:end)];
    else % no corner points in between
        ifAlreadyTeated(pos) = 1;
    end
    pos = find(ifAlreadyTeated == 0, 1);
end
result = corner_points;