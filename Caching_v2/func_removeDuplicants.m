function noDup_breakups = func_removeDuplicants(duplicate_breakups, t, demandType, leader)
%  input: duplicate_breakups 1x4 char array
%     '{C13 B15}'    
%     '{C13 B16}'   
%     '{C14 B15}'   
%     '{C14 B16}'
% output: (any three would be fine), there are 4 symbols, however the rank is 3
%     '{C13 B15}'    
%     '{C13 B16}'   
%     '{C14 B15}'  

% step 0. if duplicate_breakups is a column vector, make it a row
if size(duplicate_breakups, 2) > 1
    duplicate_breakups = duplicate_breakups';
end

% step 1. find how many unique symbols
files = demandType(leader);
posofFiles = [];
temp = duplicate_breakups{1};  % get char
for iFile = 1:length(files)
    posofFiles = [posofFiles, find(temp(1,:) == char(files(iFile)))];
end
posofFiles = sort(posofFiles);

% step 2. construct the matrix in the base of these unique symbols
numofRows = size(duplicate_breakups, 1);
symbols = [];
for irow = 1:numofRows
    temp = duplicate_breakups{irow};
    for ipos = 1:length(posofFiles)
        pos = posofFiles(ipos);
        symbols = [symbols; temp(pos:pos+t)];
    end
end
unique_symbols = unique(string(symbols)); % find the unique symbols
numofUniSymbols = length(unique_symbols);
unique_symbols = char(unique_symbols);
coefficient_matrix = zeros(numofRows, numofUniSymbols); % square matrix
for irow = 1: numofRows
    temp = duplicate_breakups{irow};
    temp = string(temp);
    for isymbol = 1:numofUniSymbols
        if ~isempty(strfind(temp, string(unique_symbols(isymbol,:))))  % e.g. find("{C14 B16}" ==  "{C14"), must be string "", can not be char ''
            coefficient_matrix(irow, isymbol) = 1;
        end
    end
end

% step 3. count rank and delete if rank < #lines (there are lines can be constructed by linear comb of others)
rankofmatrix = gfrank(coefficient_matrix, 2);
if rankofmatrix < numofRows
    noDup_breakups = duplicate_breakups(1:rankofmatrix); % keep the first (rank # of) rows
else
    noDup_breakups = duplicate_breakups;
end
