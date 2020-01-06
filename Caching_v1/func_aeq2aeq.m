function Aeq2 = func_aeq2aeq(Aeq, x)
Aeq2 = [];
begin_x = 0;
end_x = 0;
for iA = 1:length(Aeq)
    if Aeq(iA) == 0
        begin_x = begin_x + 1;
        end_x = end_x + 1;
    else % Aeq(iA) >= 1
        begin_x = begin_x + 1;
        end_x = begin_x + Aeq(iA) - 1;
        Aeq2 = [Aeq2; zeros(1, begin_x -1), ones(1, Aeq(iA)), zeros(1, length(x)-end_x)];
        begin_x = end_x; 
    end
end