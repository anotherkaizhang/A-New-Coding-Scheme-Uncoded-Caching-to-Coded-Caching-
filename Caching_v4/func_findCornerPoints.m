function corner_points = func_findCornerPoints(M_sample, R_sample)
%% find conor points
second_order_diff = diff(diff(R_sample));
corner_pos = abs(second_order_diff) > 1e-6;
corner_points = [];
ic = 1;
while ic <= length(corner_pos)
% for ic = 1:length(corner_pos)
    if corner_pos(ic) == 1
        if corner_pos(ic + 1) == 0 % pattern ... 0 1 0 ...
            corner_points = [corner_points; [M_sample(ic + 1), R_sample(ic + 1)]];
        else % pattern ... 0 1 1 0 ...
            corner_left = [M_sample(ic + 1), R_sample(ic + 1)];
            corner_left_left = [M_sample(ic), R_sample(ic)];
            corner_right = [M_sample(ic + 2), R_sample(ic + 2)];
            corner_right_right = [M_sample(ic + 3), R_sample(ic + 3)];
            corner_points = [corner_points; func_findIntersectionofTwoLines(corner_left_left, corner_left, corner_right, corner_right_right)];
            ic = ic + 1;
        end
    end
    ic = ic + 1;
end
