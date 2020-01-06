function vector = func_groupAddition(vector, num_of_ways_to_decompose_a_t_vector)
% =======================================
% Perform vector increasing without exceeding a maximum vector
% num_of_ways_in_each_TransType = [2, 2, 5],
% vector may be [1, 1, 6],
% We refine it to be [1, 2, 1],
% Find the digit that "vector" exceeds "num_of_ways_in_each_TransType", reset it to the minimum and add 1 to the previous digit,
% recursively do this until no digit in "vector" exceeds "num_of_ways_in_each_TransType"

% index = 1; % set it to "True" first
index = find((vector > num_of_ways_to_decompose_a_t_vector) == 1); % index where "vector" exceeds "num_of_ways_in_each_TransType"
while index % True
    vector(index) = 1; 
    if index > 1
        vector(index - 1) = vector(index - 1) + 1;
    else
        break
    end
    index = find((vector > num_of_ways_to_decompose_a_t_vector) == 1); % index where "vector" exceeds "num_of_ways_in_each_TransType"
end