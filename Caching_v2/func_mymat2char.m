function result = func_mymat2char(mat)
temp = char(string(mat));
temp2 = [];
for icol = 1:size(mat,2)
    temp2 = [temp2, temp(:,:,icol)];
end
result = temp2;                   % e.g. result = [3,4;3,4;3,4;3,4]
result = string(result);          % e.g. result = ["34";"34";"34";"34"]