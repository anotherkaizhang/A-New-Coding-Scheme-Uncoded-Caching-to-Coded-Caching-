% function func_partialDecomposition(N, K, t) 
clc; close all;
% if nargin == 3   % if the number of inputs equals 3
% sample_interval = 10; % then make the fourth value, sample_interval, equal to my default value, 51.
% end

%% 0.  Setting. 
N = 5; K = 5; t = 2; 
% directory = '~/home/kzhang/Dropbox/Work/3rd_paper_ISIT_2017/codes/Caching_v3/RunningTempFiles/'
directory = 'C:\Users\Kai Zhang\Documents\MATLAB\';
% directory = pwd;
files = 1 : N;

%% 1.  Plot Tian v.s. Yu.
N_tilde = min(N, K);
figure
linewidth = 1;
% Tian
t_axis = 0:K; % 
M_Tian = t_axis.*((N-1)*t_axis+K-N)/K/(K-1);
R_Tian =  N*(K - t_axis)/K;
plot(M_Tian, R_Tian, 'kd--','LineWidth',linewidth);

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
plot(M_Yu,R_Yu, 'rx-','LineWidth',linewidth);

M_Jesus = zeros(1,N);
R_Jesus = zeros(1,N);
for g = 1:N
    M_Jesus(g) = N/(K*g);
    R_Jesus(g) = N - N/K*(N+1)/(g+1);
end
plot(M_Jesus, R_Jesus, 'b+-.');

%% 2.  Demand types
partition = intpartgen(K, N);
partition = cell2mat(partition(end));

d_num = size(partition, 1); % Number of demand types

%% 4. For each demand.
num_var_total = 0;
test = 0;
for ip = 1 : size(partition, 1)
    disp([num2str(ip/size(partition, 1)*100), '%'])
    
    % Generate coefficients a: matrix "demand number" x "decomposition number in each demand"
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
    end
    t_num = size(t_vector, 1); % number of transmission types

    % All possible decompositions of a transmission type
    num_of_t_vector = size(t_vector, 1);
    num_of_ways_to_decompose_a_t_vector = zeros(num_of_t_vector, 1);
    P_d = cell(num_of_t_vector, 1);
    for i_t_vector = 1: num_of_t_vector
        present_files = find(t_vector(i_t_vector, :) > 0); % those files present at this transmission type t
        P_d{i_t_vector} = SetPartition(present_files, length(present_files)); % all possible ways to decompose this transmission type t
        num_of_ways_to_decompose_a_t_vector(i_t_vector) = length(P_d{i_t_vector});
    end

    % Choices of decomposition pattern for all transmission types
    func_generateDecompositionCombination(num_of_ways_to_decompose_a_t_vector, num_of_t_vector);

    % Construct a, R_d_pd and DeltaM_d_Pd
    num_of_total_decomposition_ways = prod(num_of_ways_to_decompose_a_t_vector);
    a_d_Pd = zeros(num_of_total_decomposition_ways + 1, 1);
    R_d_Pd = zeros(1, num_of_total_decomposition_ways + 1);
    M_d_Pd = zeros(K, num_of_total_decomposition_ways + 1);
    fileID2 = fopen('pattern_choice.txt','r');
    for ipp = 1 : num_of_total_decomposition_ways
        %% ===== This part can be deleted. The purpose is to look at the result of the LP ===
        % ====== then find the non-zero alpha, see what the corresponding partition type it is  =====
%         test = test + 1
%         if ismember(test, [1           2           4          10            24          26        1155 ...
%                            1245        1795        3117        5492 ...
%                            5493       45067       45111 ...
%                            45114       45158       45338       45435       45438       45440])
%             fprintf("demand is: ");
%             m_vector
%             fprintf("transmission types are:")
%             t_vector
%             P_d
%             keyboard
%         end
        % ============================================================================
        
        
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
    current_DeltaM = min(K-t, N_tilde) * func_mynchoosek(K-1, t-1) * ones(K, 1);
    M_d_Pd(:, end) = N * func_mynchoosek(K-1, t-1) - current_DeltaM;
    R_d_Pd(end) = min(K-t, N_tilde) * func_mynchoosek(K,t);
    
    %% ===== This part can be deleted. The purpose is to look at the result of the LP ===
    % ====== then find the non-zero alpha, see what the corresponding partition type it is  =====
%     test = test + 1
%     if ismember(test, [1           2           4          10          24          26 ...
%                      1155        1245        1795        3117        5492        5493...
%                       45067       45111       45114       45158       45338       45435...
%                        45438       45440])
%         fprintf("demand is: ");
%         m_vector
%         fprintf("special uncoded\n");
%         keyboard
%     end
    
    % ============================================================================
    
    
    % Save M_d_Pd and R_d_Pd
    num_var = length(a_d_Pd);
    filename = strcat(directory, '\N=', int2str(N),',K=',int2str(K),',t=', int2str(t), ', demand=',int2str(d), '.mat'); 
    save(filename,'M_d_Pd','R_d_Pd','P_d','a_d_Pd', 'num_var'); % ,'pattern_choice');
    num_var_total = num_var_total + num_var;
end

%% Find conor points
corner_points = [M_Tian(t+1), R_Tian(t+1); M_Yu(t+1), R_Yu(t+1)];
ifAlreadyTeated = [0]; % this segment between Tian and Yu has not been tested yet

A = zeros(K * size(partition, 1) + size(partition, 1), num_var_total);
col_current_pos = 0; % set the current position of column to 1
Aeq = zeros(size(partition, 1), num_var_total + 2);
for ip = 1 : size(partition, 1)
    disp([num2str(ip/size(partition, 1)*100), '%'])
    
    % Generate coefficients a: matrix "demand number" x "decomposition number in each demand"
    d = zeros(1, sum(partition(ip,:)));
    d_start = 1;
    for sym = 1:N  % d -> demandFile: [3,1,0](N-dimensional) -> [1,1,1,2](K-dimensional), i.e. (A,A,A,B)
        d(d_start: d_start + partition(ip, sym) -1) = files(sym)*ones(1, partition(ip, sym));
        d_start = d_start + partition(ip, sym); % the third demand type d = (A,A,B,C)
    end
    filename = strcat(directory, '\N=', int2str(N),',K=',int2str(K),',t=', int2str(t), ', demand=',int2str(d), '.mat'); 
    load(filename,'M_d_Pd','R_d_Pd','P_d','a_d_Pd', 'num_var'); % ,'pattern_choice');
    A(1 + (ip-1)*K: K + (ip-1)*K, col_current_pos + 1:col_current_pos + num_var) = M_d_Pd;  % A is a big block-diagonal matrix, each block correspond to a d
    A(K * size(partition, 1) + ip, col_current_pos + 1:col_current_pos + num_var) = R_d_Pd;
    Aeq(ip, col_current_pos + 1:col_current_pos + num_var) = ones(1, num_var);
    col_current_pos = col_current_pos + num_var;
end
A = [A, zeros(size(A, 1), 2)]; % add for M and R
A(1:K*size(partition, 1), end-1) = -nchoosek(K,t);
A(K*size(partition, 1)+1:end, end) = -nchoosek(K,t);
a_d_Pd = zeros(num_var_total + 2, 1); % reset a_d_Pd to the total number of variables, plus M and R
% f = [zeros(num_var_total, 1); (corner_points(2,2) - corner_points(1,2))/(corner_points(1,1) - corner_points(2,1)); 1];
b = zeros(K* size(partition, 1) + size(partition, 1), 1); 
beq = ones(size(partition, 1), 1);
% lb = [zeros(num_var_total, 1); corner_points(1,1); corner_points(2,2)];
% ub = [ones(num_var_total, 1); corner_points(2,1); corner_points(1,2)];
options = optimoptions('linprog','Algorithm','dual-simplex','Display','none','OptimalityTolerance',1.0000e-07);
corner_points = func_findCornerPoint(N, num_var_total, A,b,Aeq,beq, options, corner_points, ifAlreadyTeated);

legendCell = [{'Tian-Chen scheme'}; {'Yu et al. scheme'};{'Gomez-Vilardebo scheme'}];

if size(corner_points, 1) == 2
    fprintf("No corner points\n")
else
	plot(corner_points(:,1), corner_points(:,2), 'b-o')
%     legendCell = [legendCell; num2str([corner_points(2:end-1,1), corner_points(2:end-1,2)])];
    legendCell = [legendCell; {'Partial decomposition'}];
end
plot([M_Tian(t+1), M_Yu(t+1)], [R_Tian(t+1), R_Yu(t+1)], 'm-.');
legendCell = [legendCell; {'Lower convex hull of known schemes'}];
legend(legendCell, 'AutoUpdate','off')

ah = gca;
% location of the plot to be zoomed in
s_pos = [1.15,  0.65, 1.55, 1.05];
% location of the zoom-in plot
t_pos = [1.8, 0.6, 3.4, 2.1];
h=line([s_pos(1),t_pos(1)],[s_pos(4),t_pos(4)]);
set(h,'LineStyle',':','color','k');
h=line([s_pos(3),t_pos(3)],[s_pos(2),t_pos(2)]);
set(h,'LineStyle',':','color','k');
zoomPlot(ah, s_pos, t_pos);

    







