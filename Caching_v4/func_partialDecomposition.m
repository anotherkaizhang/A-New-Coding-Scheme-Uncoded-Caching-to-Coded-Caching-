% function func_partialDecomposition(N, K, t) 
clc; close all;
% if nargin == 3   % if the number of inputs equals 3
sample_interval = 10; % then make the fourth value, sample_interval, equal to my default value, 51.
% end

%% 0.  Setting. 
N = 3; K = 4; t = 2; 
% directory = '~/home/kzhang/Dropbox/Work/3rd_paper_ISIT_2017/codes/Caching_v3/RunningTempFiles/'
directory = 'C:\Users\Kai Zhang\Documents\MATLAB\';
% directory = pwd;
files = 1:N;

%% 1.  Plot Tian v.s. Yu.
N_tilde = min(N, K);
figure
% Tian
t_axis = 0:K; % 
M_Tian = t_axis.*((N-1)*t_axis+K-N)/K/(K-1);
R_Tian =  N*(K - t_axis)/K;
plot(M_Tian, R_Tian, 'k--d', 'Linewidth', 1)

% Yu 
hold on
t_axis = 0:K-1;
M_Yu = t_axis*N/K;
R_Yu = zeros(size(t_axis)); 
for iR = 1:length(t_axis)
    if K - N_tilde >= t_axis(iR)+1
        R_Yu(iR) = (nchoosek(K, t_axis(iR)+1) - nchoosek(K - N_tilde, t_axis(iR)+1)) / nchoosek(K,t_axis(iR));
    else
        R_Yu(iR) = nchoosek(K, t_axis(iR)+1) / nchoosek(K, t_axis(iR));
    end
end
M_Yu = [M_Yu, N];
R_Yu = [R_Yu, 0];
plot(M_Yu,R_Yu, 'r-x', 'Linewidth', 1)
title(['N=', num2str(N), ', K=', num2str(K)])

figure; hold on
plot(M_Tian, R_Tian, 'k--d', 'Linewidth', 1)
plot(M_Yu,R_Yu, 'r-x', 'Linewidth', 1)

%% 2.  Demand types.
partition = intpartgen(K, N);
partition = cell2mat(partition(end));
d_num = size(partition, 1); % Number of demand types

%% 3. Transmission types:
x_min = M_Tian(t+1); x_max = M_Yu(t+1);
y_min = R_Yu(t+1); y_max = R_Tian(t+1);
axis([x_min x_max y_min y_max])

% x = [0.13 0.905] y = [0.11 0.925]
x_scope = (x_max - x_min)/(0.905 - 0.13);
y_scope = (y_max - y_min)/(0.925 - 0.11);

%% 4. For each demand.
% Fix M_sample for linear programming, and compute R_sample
M_sample = linspace(M_Tian(t+1), M_Yu(t+1), sample_interval);
R_sample = zeros(size(partition, 1), sample_interval);

% profile on

for ip = 1 : size(partition, 1)
    disp([num2str(ip/size(partition, 1)*100), '%'])
    
    % Generate coefficients a: matrix "demand number" x "decomposition number in each demand"
%     a_d_Pd = [];
%     R_d_Pd = [];
%     M_d_Pd = [];

    d = zeros(1, sum(partition(ip,:)));
    d_start = 1;
    for sym = 1:N  % d -> demandFile: [3,1,0](N-dimensional) -> [1,1,1,2](K-dimensional), i.e. (A,A,A,B)
        d(d_start: d_start + partition(ip, sym) -1) = files(sym)*ones(1, partition(ip, sym));
        d_start = d_start + partition(ip, sym); % the third demand type d = (A,A,B,C)
    end
    m_vector = partition(ip, :);

    % For this t: this d: All possible transmission types
    t_vector = [];
    t_vector_min = zeros(1, N);
    t_vector_max = zeros(1, N);
    for iN = 1:N
        t_vector_min(iN) = max(0, t+1-(K-partition(ip, iN)));
        t_vector_max(iN) = min(t+1, partition(ip, iN));
    end  
    
    t_vector_candidate_tmp = intpartgen(t+1, N);
    if t+1 >= N
        t_vector_candidate = cell2mat(t_vector_candidate_tmp(end));
    else %i.e. t+1 < N
        t_vector_candidate = [cell2mat(t_vector_candidate_tmp(end)), zeros(size(cell2mat(t_vector_candidate_tmp(end)), 1), N - (t+1))];
    end
    
    for icand = 1:size(t_vector_candidate, 1)
        temp = perms(t_vector_candidate(icand, :)); % for each candidate, generate all permutations
        temp = unique(temp, 'rows'); % might exist repetation
        satisfy_cond = (sum(t_vector_min <= temp, 2) + sum(temp <= t_vector_max, 2) == 2*N);
        t_vector_tmp = temp(find(satisfy_cond), :);
        t_vector = [t_vector; t_vector_tmp];
%         for itemp = 1:size(temp, 1)
%             if sum(t_vector_min <= temp(itemp,:)) + sum(temp(itemp,:) <= t_vector_max) == 2*N
%                 t_vector = [t_vector; temp(itemp, :)];
%             end
%         end

    end
    t_num = size(t_vector, 1); % number of transmission types

    % All possible decompositions of a transmission type
    num_of_t_vector = size(t_vector, 1);
    num_of_ways_to_decompose_a_t_vector = zeros(num_of_t_vector, 1);
    P_d = cell(num_of_t_vector, 1);
    for i_t_vector = 1: num_of_t_vector
        present_files = find(t_vector(i_t_vector, :) > 0); % those files present at this transmission type t
        P_d{i_t_vector} = SetPartition(present_files); % all possible ways to decompose this transmission type t
        num_of_ways_to_decompose_a_t_vector(i_t_vector) = length(P_d{i_t_vector});
    end

    % Choices of decomposition pattern for all transmission types
    func_generateDecompositionCombination(num_of_ways_to_decompose_a_t_vector, num_of_t_vector);

    % Construct a, R_d_pd and DeltaM_d_Pd
    num_of_total_decomposition_ways = prod(num_of_ways_to_decompose_a_t_vector);
    a_d_Pd = zeros(1, num_of_total_decomposition_ways + 1);
    R_d_Pd = zeros(1, num_of_total_decomposition_ways + 1);
    M_d_Pd = zeros(K, num_of_total_decomposition_ways + 1);
    fileID2 = fopen('pattern_choice.txt','r');
    for ipp = 1 : num_of_total_decomposition_ways
        % pattern = pattern_choice(:, ipp);  % fix a partition pattern for all t_vector
        pattern = fscanf(fileID2, '%d', size(num_of_ways_to_decompose_a_t_vector));
        % Calculate current_R
        current_R = 0;
        for t_vector_index = 1 : t_num % for each t_vector
            current_t_vector = t_vector(t_vector_index, :);
            P_t_d = P_d{t_vector_index}{pattern(t_vector_index)}; % all partitions in a t_vector
            for iP = 1 : size(P_t_d, 2) % how many partitions in this t_vector
                P = P_t_d{iP};
                temp1 = 1; temp2 = 1; temp3 = 1;
                for in = 1: length(P) % for each file in P
                    n = P(in);
                    temp1 = temp1 * func_mynchoosek(m_vector(n), current_t_vector(n));
                    temp2 = temp2 * func_mynchoosek(m_vector(n)-1, current_t_vector(n));
                end
                supp_t = find(current_t_vector > 0);
                suppt_exclude_P = setdiff(supp_t, P);
                if ~isempty(suppt_exclude_P)
                    for is = 1 : length(suppt_exclude_P)
                        n = suppt_exclude_P(is);
                        temp3 = temp3 * func_mynchoosek(m_vector(n),  current_t_vector(n));
                    end
                end
                current_R = current_R + (temp1 - temp2) * temp3;
            end
        end

        % Calculate current_DeltaM
        current_DeltaM = zeros(K, 1);
        for k = 1:K
            current_DeltaM_k = 0;
            for t_vector_index = 1 : t_num % for each t_vector
                current_t_vector = t_vector(t_vector_index, :);
                if current_t_vector(d(k)) <= 0 % condition: t_(d_k) > 0
                    continue;
                end
                P_t_d = P_d{t_vector_index}{pattern(t_vector_index)}; % all partitions in a t_vector
                temp1 = func_mynchoosek(m_vector(d(k)) - 1, current_t_vector(d(k)) - 1);
                current_DeltaM_k_inside = 0;
                for iP = 1 : size(P_t_d, 2) % how many partitions in this t_vector
                    P = P_t_d{iP};
                    if ~isempty(find(P == d(k), 1))  % condition: d_k \notin P
                        continue;
                    end
                    temp2 = 1; temp3 = 1; temp4 = 1;
                    for in = 1:length(P)
                        n = P(in);
                        temp2 = temp2 * func_mynchoosek(m_vector(n), current_t_vector(n));
                        temp3 = temp3 * func_mynchoosek(m_vector(n)-1, current_t_vector(n));
                    end
                    supp_t = find(current_t_vector > 0);
                    suppt_exclude_P_and_dk = setdiff(supp_t, union(P, d(k)));
                    if ~isempty(suppt_exclude_P_and_dk)
                        for is = 1 : length(suppt_exclude_P_and_dk)
                            n = suppt_exclude_P_and_dk(is);
                            temp4 = temp4 * func_mynchoosek(m_vector(n), current_t_vector(n));
                        end
                    end
                    current_DeltaM_k_inside = current_DeltaM_k_inside + (temp2 - temp3) * temp4;
                end
                current_DeltaM_k = current_DeltaM_k + temp1 * current_DeltaM_k_inside;
            end
            current_DeltaM(k) = current_DeltaM_k;
        end

        R_d_Pd(ipp) = current_R;
        M_d_Pd(:,ipp) = N*nchoosek(K-1,t-1) - current_DeltaM;
    end
    fclose(fileID2);
    % Add special uncoded transmission
%     if K - t >= N_tilde
%         a_d_Pd = [a_d_Pd, 0];
        current_DeltaM = min(K-t, N_tilde) * func_mynchoosek(K-1, t-1) * ones(K, 1);
        M_d_Pd(:, end) = N * func_mynchoosek(K-1, t-1) - current_DeltaM;
        R_d_Pd(end) = min(K-t, N_tilde) * func_mynchoosek(K,t);
%     end

    % Save M_d_Pd and R_d_Pd
    filename = strcat(directory, '\N=', int2str(N),',K=',int2str(K),',t=', int2str(t), ', demand=',int2str(d), '.mat'); 
    save(filename,'M_d_Pd','R_d_Pd','P_d','a_d_Pd'); % ,'pattern_choice');
    
    % Linear programming
    num_var = length(a_d_Pd);
    for im = 1:sample_interval
        f = R_d_Pd';
        A = M_d_Pd;
        b = M_sample(im) * nchoosek(K,t) * ones(K, 1);
        Aeq = ones(1, num_var);
        beq = 1;
        lb = zeros(num_var, 1);
        ub = ones(num_var, 1);
        options = optimoptions('linprog','Algorithm','dual-simplex','Display','none','OptimalityTolerance',1.0000e-07);
        [a_d_Pd, fval] = linprog(f,A,b,Aeq,beq,lb, ub, options);
        R_sample(ip, im) = a_d_Pd' * R_d_Pd' / nchoosek(K,t);
    end
end
% profile off
% profile viewer

%% Add function (Dr.Tian's request): For full demands, label all corner points
% full_demand = find(partition(:, end) > 0);
% for ifu = full_demand' % full_demand'
%     d = [];
%     for sym = 1:N  % d -> demandFile: [3,1,0](N-dimensional) -> [1,1,1,2](K-dimensional), i.e. (A,A,A,B)
%         d = [d, files(sym)*ones(1, partition(ifu, sym))]; % the third demand type d = (A,A,B,C)
%     end
%     corner_points = func_findCornerPoints(M_sample, R_sample(ifu, :));
%     for ic = 1:size(corner_points, 1)
%         filename = strcat(directory, '\N=', int2str(N),',K=',int2str(K),',t=', int2str(t), ', demand=', int2str(d), '.mat'); 
%         load(filename,'M_d_Pd','R_d_Pd','P_d','a_d_Pd'); % ,'pattern_choice');
%         num_var = length(a_d_Pd);
%         f = R_d_Pd';
%         A = M_d_Pd;
%         b = corner_points(ic, 1) * nchoosek(K, t) * ones(K, 1);
%         Aeq = ones(1, num_var);
%         beq = 1;
%         lb = zeros(num_var, 1);
%         ub = ones(num_var, 1);
%         options = optimoptions('linprog','Algorithm','dual-simplex','Display','none','OptimalityTolerance',1.0000e-07);
%         [a_d_Pd, fval] = linprog(f,A,b,Aeq,beq, lb, ub, options);
%         corner_point_partition = find(a_d_Pd > 0);

%         text = [];
%         for il = 1:length(corner_point_partition)
%             partition_index_in_each_t = pattern_choice(:, corner_point_partition(il));
%             text1 = []; % text to be print aside the corner point
%             for l = 1:size(partition_index_in_each_t, 1) % decomposition pattern in each t
%                 cellarray = P_d{l}{partition_index_in_each_t(l)};
%                 text2 = [];
%                 for icell = 1:length(cellarray)
%                     if length(cell2mat(cellarray(icell))) == 1
%                         text2 = strcat(text2, '[', mat2str(cell2mat(cellarray(icell))), ']');
%                     else
%                         text2 = strcat(text2, mat2str(cell2mat(cellarray(icell))));
%                     end
%                     if icell < length(cellarray)
%                         text2 = strcat(text2, ','); % put ',' to separate each partition
%                     end
%                 end
%                 text1 = vertcat(text1, {text2}); % put '#' to separate each t
%             end
%             text = vertcat(text, {'------'});
%             text = vertcat(text, text1);
%         end
%         text = sprintf('%s\n',text{:});
% 
%         x = [0.10+(corner_points(ic, 1) - x_min)/x_scope, 0.13+(corner_points(ic, 1) - x_min)/x_scope];
%         y = [0.08+(corner_points(ic, 2) - y_min)/y_scope, 0.11+(corner_points(ic, 2) - y_min)/y_scope];
%         annotation('textarrow',x,y,'String',text)
%     end
% end

plot(M_sample, R_sample)
R_sample = max(R_sample);
plot(M_sample, R_sample, '--r', 'LineWidth', 2)
legendCell = [{'Tian'}; {'Yu'}];
legendCell = [legendCell; cellstr(strcat(repmat('demand = ', d_num, 1), num2str(partition, '(%d %d %d)')))];
legendCell = [legendCell; 'Intersection'];

%% Find conor points
corner_points = func_findCornerPoints(M_sample, R_sample);
if isempty(corner_points)
    fprintf("No corner points\n")
else
	plot(corner_points(:,1), corner_points(:,2), 'r*')
    legendCell = [legendCell; num2str([corner_points(:,1), corner_points(:,2)])];
end
legend(legendCell)
title(['N=', num2str(N), ', K=', num2str(K), ', t=', num2str(t)])





