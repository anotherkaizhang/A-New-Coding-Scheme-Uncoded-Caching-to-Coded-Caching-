function dispPartObj = func_waysToBreakUp(aLinearComb)
    iaL = 2;
    while iaL <= length(aLinearComb) 
        for jaL = 1:iaL-1
            if strncmpi(aLinearComb{iaL}, aLinearComb{jaL}, 1)
                aLinearComb{jaL} = strcat(aLinearComb{jaL}, aLinearComb{iaL});
                aLinearComb{iaL} = [];
                break;
            end
        end
        iaL = iaL + 1;
    end
    empties = find(cellfun('isempty',aLinearComb)); % identify the empty cells
    aLinearComb(empties) = [];                      % remove the empty cells
    p = SetPartition(length(aLinearComb));
    dispPartObj = DispPartObj(p, aLinearComb);