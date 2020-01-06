function [ub, beq] = func_calcUbBeq(startPointx, endPointx, Aeq, beq)
ub = endPointx - startPointx;
if ~isempty(Aeq)
    for iA = 1:size(Aeq,1)
        posof1 = find(Aeq(iA,:) == 1);
        ub(posof1) = sum(endPointx(posof1) - startPointx(posof1)); % each one of the four should not exceed a number
        beq(iA) = sum(endPointx(posof1) - startPointx(posof1)); % sum change of the four should not exceed a number
    end
end
% this is to force the negative value (e.g. -1.6434E-16) to be zero
ub(ub<0)=0; 
