function [num_variables, M, R] = func_createX(groupPartationArray, t, demandType, leader)
%  input: groupPartationArray
%  e.g.   7×4 cell array
%     '{C13 B15} {A35}'    '{C13 A35} {B15}'    '{C13} {B15 A35}'    '{C13} {B15} {A35}'  
%     '{C13 B16} {A36}'    '{C13 A36} {B16}'    '{C13} {B16 A36}'    '{C13} {B16} {A36}'  
%     '{C14 B15} {A45}'    '{C14 A45} {B15}'    '{C14} {B15 A45}'    '{C14} {B15} {A45}'  
%     '{C14 B16} {A46}'    '{C14 A46} {B16}'    '{C14} {B16 A46}'    '{C14} {B16} {A46}' 
%     '{C23 B25} {A35}'    '{C23 A35} {B25}'    '{C23} {B25 A35}'    '{C23} {B25} {A35}'  
%     '{C23 B26} {A36}'    '{C23 A36} {B26}'    '{C23} {B26 A36}'    '{C23} {B26} {A36}'  
%     '{C24 B25} {A45}'    '{C24 A45} {B25}'    '{C24} {B25 A45}'    '{C24} {B25} {A45}'  
%           x_1                    x_2                 x_3                   x_4          x_1 + x_2 + x_3 + x_4 = 1
%  result:
%     num_variabes = 4;
%     M =                                [x_1]
%  user1 [                       ]    x =[x_2]  R = [  ,  ,  ,  ]
%  user2 [                       ]       [x_3]
%  user3 [                       ]       [x_4]
%  user4 [                       ]    
%  user5 [                       ]
%  user6 [                       ]

num_variables = size(groupPartationArray, 2);        
for x_i = 1:num_variables % for each x_i
%     haha = cell(1,3);
%     haha{1,1} = '{C13 B15} {A35}'
%     haha{1,2} = '{C13 B16} {A36}'
%     haha{1,3} = '{C15 B16} {A36}'
%     cells = func_seperateFiles(haha, t, demandType, leader); % for x_1
    cells = func_seperateFiles(groupPartationArray(:, x_i), t, demandType, leader); % for x_1
    [M(:,x_i), R(x_i)] = func_countRank(cells, demandType, leader);
    R(x_i) = R(x_i) - size(groupPartationArray,1);
end
M = M'; % each user is a column