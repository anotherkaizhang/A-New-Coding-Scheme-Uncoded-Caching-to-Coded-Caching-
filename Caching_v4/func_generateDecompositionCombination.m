% function pattern_choice = func_generateDecompositionCombination(num_of_ways_to_decompose_a_t_vector, num_of_t_vector)
function func_generateDecompositionCombination(num_of_ways_to_decompose_a_t_vector, num_of_t_vector)
% ==============================================
% e.g. num_of_ways_in_each_TransType = [2, 2, 5]
%      means there are 2 ways to decompose transmission type 1,
% .                    2 ways to decompose transmission type 2,
%                      5 ways to decompose transmission type 3,
% .    there are 2*2*5=20 ways in total, 
% .    the result should be: [1, 1, 1;
%                            [1, 1, 2;
%                               ....
%                            [2, 2, 5]


% pattern_choice = zeros(num_of_t_vector, prod(num_of_ways_to_decompose_a_t_vector)); % this matrix is too big, will generate error when running (5,11) case
fileID = fopen('pattern_choice.txt','w');
start = ones(num_of_t_vector, 1); % start from [1, 1, 1]
fprintf(fileID,'%d\r\n',start);
% pattern_choice(:, 1) = start; % set the first one to [1, 1, 1]
ind = 2;
addition = zeros(num_of_t_vector, 1); addition(end) = 1;
% lastDigit = num_of_TransType; % last digit 
while ~isequal(start, num_of_ways_to_decompose_a_t_vector)
    start = func_groupAddition(start + addition, num_of_ways_to_decompose_a_t_vector);
    fprintf(fileID,'%d\r\n',start);
    % pattern_choice(:, ind) = start;
    ind = ind + 1;
end
fclose(fileID);