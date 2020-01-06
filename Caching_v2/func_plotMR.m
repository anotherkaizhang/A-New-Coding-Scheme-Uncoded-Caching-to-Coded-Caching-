function func_plotMR(P_all, numofPartition)  
%% plot new (M,R) points
opera = 'intersection';
newPointPartation = [];
% find which demand generates new points
for ip = 1:numofPartition 
    if size(P_all{ip}, 1) >= 3
        newPointPartation = [newPointPartation, ip];
    end
end

if length(newPointPartation) == 0 % no demand generates new points
    % do nothing 
elseif length(newPointPartation) == 1
    temp = P_all{newPointPartation}; temp = [temp; temp(1,:)]; xb = temp(:,1); yb = temp(:,2);
    xb(end) = []; yb(end) = [];
    plot(xb, yb, 'g-d', 'Linewidth', 2) 
else
    temp = P_all{newPointPartation(1)}; temp = [temp; temp(1,:)]; xb = temp(:,1); yb = temp(:,2);
    for inew = 2:length(newPointPartation)
        temp = P_all{newPointPartation(inew)}; temp = [temp; temp(1,:)]; x1 = temp(:,1); y1 = temp(:,2);
        [xb, yb] = polybool(opera, xb, yb, x1, y1);
    end
    xb(end) = []; yb(end) = [];
    plot(xb, yb, 'g-d', 'Linewidth', 2) 
end

    
    