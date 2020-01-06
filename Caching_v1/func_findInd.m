function x = func_findInd(x, Aeq2)
% =====================================
% find the indices of x that should be deleted,
% e.g. 
%      x =  1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1
%   Aeq2 = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0;
%           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]
% then x = 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1 
% all indices are deleted (reaches Tian point)
% =====================================
if isempty(Aeq2)
    ind = find(x == 1);
else
    for irow = 1:size(Aeq2, 1)
        posof1 = find(Aeq2(irow,:) == 1);
        if ~isempty(find(x(posof1) == 1))
           x(posof1) = 1; 
        end
    end
end


