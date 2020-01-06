function x_pos = func_xPos(dispPartObj, file, t)
x_pos = [];
fileApos = find(dispPartObj{1} == file); % where 'A' appears
firstPartition = dispPartObj{1}; % '{B234B235 A245A345}'
filegroup = firstPartition(fileApos(1): fileApos(end) + t); % A245A345
filegroup = strcat('{',filegroup,'}') % {A245A345}
for id = 1:length(dispPartObj)
    if strfind(dispPartObj{id}, filegroup);
        x_pos = [x_pos, id];
    end
end
x_pos = x_pos - 1;