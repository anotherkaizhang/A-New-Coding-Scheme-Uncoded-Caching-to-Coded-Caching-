function group = func_findGroups(linearCombFile, t)
% 
%   linearCombFile: 19×3 string array
% 
%     "66"    "65"    "65"
%     "66"    "65"    "65"
%     "67"    "65"    "65"
%     "67"    "65"    "65"
%     "66"    "66"    "65"
%     "67"    "66"    "65"
%     "67"    "66"    "65"
%     "67"    "66"    "65"
%     "67"    "66"    "65"
%     "67"    "67"    "65"
%     "66"    "66"    "65"
%     "67"    "66"    "65"
%     "67"    "66"    "65"
%     "67"    "66"    "65"
%     "67"    "67"    "65"
%     "67"    "66"    "66"
%     "67"    "66"    "66"
%     "67"    "67"    "66"
%     "67"    "67"    "66"
% groups: 1x7 cell array
%   {[1;2],[5;11],[3;4],[6;7;8;9;12;13;14],[16;17],[10;15],[18;19]}

a = char(linearCombFile);
temp = [];
for it = 1:t+1
    temp = [temp, a(:,:,it)];
end
temp = str2num(temp);
unique_sets = unique(temp);
group = {}
for iu = 1:length(unique_sets)
    % delete those can not break up, e.g. A12+A13+A23
    ind = find(temp == unique_sets(iu));
    if length(unique(linearCombFile(ind(1),:))) > 1
        group{iu} = find(temp == unique_sets(iu));
    end
end
group = group(~cellfun('isempty',group));
