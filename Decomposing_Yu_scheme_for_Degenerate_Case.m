clc ;close all; clear all

%% All files are requested:
fprintf('All files are requested: \n')
N = 4; K = 8; t = 2;
mn = [3,2,1];
% mn = [4,1,1];

tn = [3,2,1,0];
sum1 = 0;
for ii = 1:length(tn)
    if (K-mn(1) >= t+1-tn(ii))&(t+1-tn(ii) >=0) && (mn(1)-1 >=tn(ii)-1 )&(tn(ii)-1 >= 0)
        temp = nchoosek(K-mn(1), t+1-tn(ii))*nchoosek(mn(1)-1,tn(ii)-1);
        sum1 = sum1 + temp;
        fprintf('The total number of transmitted symbols for tn = %d is %d \n', tn(ii), temp)
    end
end
fprintf('LHS = %d, RHS = %d, check if match %d \n', sum1,  nchoosek(K-1,t), sum1==nchoosek(K-1,t))

sum2 = 0;
for ii = 2:N
    for jj = 1:length(tn)
        if (K-mn(ii)-1 >= t-tn(jj))&(t-tn(jj) >=0) && (mn(ii)-1 >=tn(jj)-1 )&(tn(jj)-1 >= 0)
            temp = nchoosek(K-mn(ii)-1, t-tn(jj))*nchoosek(mn(ii)-1, tn(jj)-1);
            fprintf('The number of useful symbols collected for file %d, tn = %d is %d \n', ii, tn(jj), temp)
            sum2 = sum2 + temp;
        end
    end
end
fprintf('LHS = %d, RHS = %d, check if match %d \n', sum2,  (N-1)*nchoosek(K-2,t-1), sum2 == (N-1)*nchoosek(K-2,t-1))

%% Some files are requested:
fprintf('Only some files are requested: \n')
N = 3; K = 6; t = 2; N_star = 2;
mn = [3,3,0];
tn_star = [3,2,1,0];
sum1 = 0;
for ii = 1:length(tn_star)
    if (K-mn(1) >= t+1-tn_star(ii))&(t+1-tn_star(ii) >=0) && (mn(1)-1 >=tn_star(ii)-1 )&(tn_star(ii)-1 >= 0)
        temp = nchoosek(K-mn(1), t+1-tn_star(ii))*nchoosek(mn(1)-1,tn_star(ii)-1);
        sum1 = sum1 + temp;
        fprintf('The total number of transmitted symbols for tn = %d is %d \n', tn_star(ii), temp)
    end
end
fprintf('LHS = %d, RHS = %d, check if match %d \n', sum1,  nchoosek(K-1,t), sum1==nchoosek(K-1,t))

sum2 = 0;
for ii = 2:N_star
    for jj = 1:length(tn_star)
        if (K-mn(ii)-1 >= t-tn_star(jj))&(t-tn_star(jj) >=0) && (mn(ii)-1 >=tn_star(jj)-1 )&(tn_star(jj)-1 >= 0)
            temp = nchoosek(K-mn(ii)-1, t-tn_star(jj))*nchoosek(mn(ii)-1, tn_star(jj)-1);
            fprintf('The number of useful symbols collected for file %d, tn = %d is %d \n', ii, tn_star(jj), temp)
            sum2 = sum2 + temp;
        end
    end
end
sum2 = sum2 + (N - N_star)*nchoosek(K-2, t-1);
fprintf('LHS = %d, RHS = %d, check if match %d \n', sum2,  (N-1)*nchoosek(K-2,t-1), sum2 == (N-1)*nchoosek(K-2,t-1))
