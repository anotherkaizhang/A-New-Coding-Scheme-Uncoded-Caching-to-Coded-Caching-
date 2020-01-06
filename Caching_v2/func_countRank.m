function [M, R] = func_countRank(cells, demandType, leader)
%  input: duplicate_breakups 1x3 char array, already removes duplications
%     '{C13 B15}'      '{A35}'
%     '{C13 B16}'   OR '{A36}'
%     '{C14 B15}'      '{A45}'
%                      '{A46}'
%        acell          anothercell
%    above are two cells,
%  output:              
% M = user1 user2 user3 user4 user5 user6
%    [                                   ];
% R = [];
K = length(demandType);
M = zeros(1,K); R = 0;
for icell = 1:length(cells)
    acell = cells{icell};
    % step 0:
    if size(acell, 2) > 1
        acell = acell';
    end
    % Step 1: for each '{C13 B15}' or '{A45}', find the user: 1 OR 4,5, then add a to
    % the corresponding places in M,
    numofRows = size(acell, 1);
    files = demandType(leader);
    numofUsers = length(demandType);
    temp = acell{1}; % acell(1) gives cell, acell{1} gives char
    numofFiles = 0;
    for iFile = 1:length(files)
        numofFiles = numofFiles + sum(temp(1,:) == char(files(iFile)));
    end
    howmanynum = 0; % find how many numbers equal to files
    for iUser = 1:numofUsers
        if sum(temp(1, :) == int2str(iUser)) == numofFiles
            howmanynum = howmanynum + 1;
        end
    end
    index = zeros(numofRows, howmanynum); % find whay number equal to files                      
    for irow = 1:numofRows
        numbers = [];
        temp = acell{irow};
        for iUser = 1:numofUsers
            if sum(temp == int2str(iUser)) == numofFiles
                numbers = [numbers, iUser];
            end
        end
        index(irow,:) = numbers;
    end
    for iInd = 1:size(index, 1)
        M(index(iInd,:) )=  M(index(iInd,:)) + 1;
    end
    R = R + size(index, 1);
end


