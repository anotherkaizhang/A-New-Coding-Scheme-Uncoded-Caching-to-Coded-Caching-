function cells = func_seperateFiles(samePartationFromaGroup, t, demandType, leader)
% input: samePartationFromaGroup: cell array
%     '{C13 B15} ||  {A35}'     
%     '{C13 B16} ||  {A36}'    
%     '{C14 B15} ||  {A45}'    
%     '{C14 B16} ||  {A46}'
%     ============
%     '{C23 B25} ||  {A35}'      
%     '{C23 B26} ||  {A36}'     
%     '{C24 B25} ||  {A45}' 
%        col{1}     col{2}
% return:
%     as seperated

% step 0: if samePartationFromaGroup is a column vector, make it a row
if size(samePartationFromaGroup, 2) > 1
    samePartationFromaGroup = samePartationFromaGroup';
end
% step 1: separate vertically using ||:
% using the first  '{C13 B15} ||  {A35}' , find '{' to break up 
files = demandType(leader);
numofUsers = length(demandType);
begin_cur = find(samePartationFromaGroup{1} == '{'); % this works since '{' is just one symbol, if multiple, must use string
end_cur = find(samePartationFromaGroup{1} == '}');
numofRows = size(samePartationFromaGroup, 1);
for ib = 1:length(begin_cur) % every 'column'
    b = begin_cur(ib); e = end_cur(ib); % begin & end of this column
    temp = char(zeros(numofRows, e-b+1));
    for irow = 1:numofRows
        abreakup = samePartationFromaGroup{irow};
        temp(irow,:) = abreakup(b:e);
    end
    col{ib} = temp;
end

% till now, each col{i} is a char array

% step 2: break up each col{i},
cells = {};
c = 0;
for ib = 1:length(begin_cur) % each column
    temp = col{ib}; % CLASS: char array
    if length(temp(1,:)) == 3 + t % '{' , 'A' , '}' % the second column, just find unique ones
        c = c + 1;
        cells{c} = cellstr(char(unique(string(temp))));  % stores cell array
    else  % the first column, we have to break up by counting ranks 
        % count how many files
        numofFiles = 0;
        for iFile = 1:length(files)
            numofFiles = numofFiles + sum(temp(1,:) == char(files(iFile)));
        end
        % find what number equal to files
        howmanynum = 0;
        for iUser = 1:numofUsers
            if sum(temp(1, :) == int2str(iUser)) == numofFiles
                howmanynum = howmanynum + 1;
            end
        end
        index = zeros(numofRows, howmanynum);                            
        for irow = 1:numofRows
            numbers = [];
            for iUser = 1:numofUsers
                if sum(temp(irow, :) == int2str(iUser)) == numofFiles
                    numbers = [numbers, iUser];
                end
            end
            index(irow,:) = numbers;
        end
        % seperate according to numbers
        [~,~,label] = unique(func_mymat2char(index)); % if 
        for iL = 1:max(label)
            c = c+1;
            whichrow = find(label == iL);
%             haha = cell(1,4);
%             haha{1,1} = '{C13 B15}'   ; 
%             haha{1,2} = '{C13 B16}'  ; 
%             haha{1,3} = '{C14 B15}' ;  
%             haha{1,4} = '{C14 B16}';
%             cells{c} = func_removeDuplicants(haha, t, demandType, leader);
            cells{c} = func_removeDuplicants(cellstr(temp(whichrow, :)), t, demandType, leader); % function only takes cell array
        end
    end
end
