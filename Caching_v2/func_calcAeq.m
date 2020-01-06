function Aeq = func_calcAeq(num_variables, num_variables_total)
Aeq = [];
begin_x = 0;
end_x = 0;
for iA = 1:length(num_variables)
    if num_variables(iA) == 1
        begin_x = begin_x + 1;
        end_x = end_x + 1;
    else % Aeq(iA) >= 1
        begin_x = begin_x + 1;
        end_x = begin_x + num_variables(iA) - 1;
        Aeq = [Aeq; zeros(1, begin_x -1), ones(1, num_variables(iA)), zeros(1, num_variables_total - end_x)];
        begin_x = end_x; 
    end
end