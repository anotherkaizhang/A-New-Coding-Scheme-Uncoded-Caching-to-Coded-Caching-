function demandName = func_findBetaNew(M_matrix, R_matrix, alpha, demandNum, K)
%             d1   d2   d3
%  demand = [ M    M    M ]
%           [ R    R    R ]
%
%
demand = zeros(2, demandNum); 
for id = 1:demandNum
    demand(1, id) = min(M_matrix(1+(id-1)*K: id*K));
    demand(2, id) = alpha*R_matrix(id);
end
demand = sum(demand); % find the user letting (alpha*R + M) to be the maximum
demandName = find(demand == max(demand));
% if all 3 demands reaches maximum, pick any one
demandName = demandName(1); 

