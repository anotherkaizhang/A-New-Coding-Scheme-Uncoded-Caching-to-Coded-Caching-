function x = func_calcEndPointx(Aeq_, x_)
x = ones(size(x_)); 
if ~isempty(Aeq_) % if there is \sum = 1 condition
    for iA = 1:size(Aeq_, 1)
        posof1 = find(Aeq_(iA,:) == 1);
        x(posof1(1:end-1)) = 0;
    end
end