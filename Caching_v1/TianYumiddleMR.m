% function plotMR(N, K, Ned)
clc; close all; clear all

%% 0. Setting: 
N = 3; K = 5;
A = sym('A'); 
B = sym('B');
C = sym('C');
D = sym('D');
E = sym('E');
F = sym('F');
G = sym('G');
H = sym('H');
I = sym('I');
J = sym('J'); 
% K = sym('K')
L = sym('L');
M = sym('M');
% N = sym('N')
O = sym('O');
P = sym('P');
Q = sym('Q');
R = sym('R');
S = sym('S');
T = sym('T');
U = sym('U');
V = sym('V');
W = sym('W');
X = sym('X');
Y = sym('Y');
Z = sym('Z');
files = [A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z];

% 0.1.  Plot Tian v.s. Yu:
Ned = min(K,N); % only consider demand types where all files are requested
figure
% Tian
t_axis = 0:K;
M_Tian = t_axis.*((N-1)*t_axis+K-N)/K/(K-1);
R_Tian =  N*(K-t_axis)/K;
plot(M_Tian,R_Tian, 'm-*', 'Linewidth', 2) 

% Yu 
hold on
t_axis = 0:K-1;
M_Yu = t_axis*N/K;
R_Yu = zeros(size(t_axis)); 
for iR = 1:length(t_axis)
    if K - Ned >= t_axis(iR)+1
        R_Yu(iR) = (nchoosek(K, t_axis(iR)+1) - nchoosek(K - Ned, t_axis(iR)+1)) / nchoosek(K,t_axis(iR));
    else
        R_Yu(iR) = nchoosek(K, t_axis(iR)+1) / nchoosek(K, t_axis(iR));
    end
end
M_Yu = [M_Yu, N];
R_Yu = [R_Yu, 0];
plot(M_Yu,R_Yu,'b-o', 'Linewidth', 2) 
legend('Tian','Yu')

% 0.2 Demand types:
partition = intpartgen(K - N);
partition = cell2mat(partition(end));
if size(partition, 2) < N
    partition = [partition, zeros(size(partition,1), N - size(partition, 2))];
end
partition = ones(size(partition)) + partition; % all files are requested: each file at least one user

%% 1. Iterate each t (=2 for testing)
t_axis = 0:K;
for it = 2:length(t_axis) - 1 
    start_point = [M_Yu(it), R_Yu(it)]; % M_Yu
    end_point = [M_Tian(it), R_Tian(it)]; % Tian
    t = t_axis(it);
    %% 2. Iterate each demand (= (A,A,B,B,C,C) for testing):
    for ip = 1:size(partition, 1)
        AA = []; F = []; x = []; Aeq = [];
        demandType = [];
        for sym = 1:N
            demandType = [demandType, files(sym)*ones(1, partition(ip, sym))]; 
        end
        leader = zeros(1,N);
        for ifile = 1:N % leader for each file
            leader(ifile)= find(demandType == files(ifile),1);
        end
        
        
        % 2.1. t+1 transmissions (label of Y):
        t_plus_1_subsets = nchoosek(1:K,t+1); % all t+1 subsets
        % 2.2. Delete redundant Y:
        if length(leader) < t+1
            fprintf('There is no redundant linear combinations')
        else
            non_leader = setxor(1:K,leader);
            licombDel = nchoosek(non_leader,t+1); % linear combinations ready to delete
            [t_plus_1_subsets,ia] = setdiff(t_plus_1_subsets, licombDel, 'rows');
        end
        
        
        % 2.3  Symbols in each Y:
        num_subsets = size(t_plus_1_subsets,1); % # of Y
        linearComb = strings([num_subsets, t+1]); % initialize linear combination
        linearComb2 = linearComb; % This is a copy
        linearComb3 = linearComb; % This is another copy
        for iLiComb = 1:num_subsets % each Y
            t_plus_1_subset = t_plus_1_subsets(iLiComb,:);
            t_subset = nchoosek(t_plus_1_subset,t); % each symbol in this Y
            % 2.3.1  Create symbols in this Y:
            for iSymbol = 1:length(t_subset) 
                filename = demandType(setdiff(t_plus_1_subset, t_subset(iSymbol,:)));
                symbol = strcat(char(filename),strjoin(string(t_subset(iSymbol,:)),''));
                linearComb(iLiComb, iSymbol) = symbol; % this is for finding groups, will be ruined
                linearComb2(iLiComb, iSymbol) = char(filename); % this is for checking singles, will be ruined
                linearComb3(iLiComb, iSymbol) = symbol; % this is the original
            end
            % 2.3.2 Keep only the file symbol appearing once in this Y:
            for ifile = 1:N
                filename = files(ifile);
                pos = find(linearComb2(iLiComb, :) == char(filename));
                if length(pos) >=2
                    linearComb(iLiComb, pos) = randn(1,length(pos)); % set to random to prevent duplicates
                end
            end
        end
        
        
        % 2.4  Ralationship matrix (for finding groups):
        [un idx_last idx] = unique(linearComb, 'stable'); % 'un': unique symbols
        unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {sort(x)}); % positions of each unique symbol
        relationship = zeros(num_subsets, num_subsets); % initialize relationship matrix
        for iu = 1:length(unique_idx)
            if length(cell2mat(unique_idx(iu))) > 1  % pick out the repeated symbols 
                tmp = mod(cell2mat(unique_idx(iu)), size(linearComb, 1));
                tmp(tmp == 0) = num_subsets; 
                unique_idx{iu} = tmp;
                tmp = nchoosek(tmp,2);
                row = tmp(:,1) ;
                col = tmp(:,2);
                relationship(sub2ind(size(relationship), row, col)) = 1;
            else
            end
        end
        
        
        % 2.5  Find the lonely groups (linear combinations), and AA, x, f:
        lonelyGroup = find(sum(relationship + relationship', 2) == 0);
        numXcreated = zeros(length(lonelyGroup), 2);
        for ig = 1:length(lonelyGroup) % each lonely group
            aLinearComb = cellstr(linearComb3(lonelyGroup(ig),:));
            % 2.5.1  In each Y, put the same file together, and find all possible ways to break up:
            dispPartObj = func_waysToBreakUp(aLinearComb);
            num_breakup = length(dispPartObj); % # of break-ups for this Y
            if num_breakup > 1 % for each break-up, create a variable x, except the first one (do no break)
                a = zeros(num_breakup - 1, K); f = zeros(1, num_breakup - 1);
                for ip = 2:num_breakup
                    breakup = dispPartObj{ip};
                    [a(ip-1,:), f(ip-1)] = func_genOptLonelyGroup(breakup, files, N, K, t);
                end
                aeq = (num_breakup > 2) * (num_breakup - 1); % more than one way to breakup, sum must be 1
                AA = [AA; a]; F = [F, f]; x = [x; zeros(num_breakup - 1, 1)]; Aeq = [Aeq, aeq];
                clear a, clear f
            end
            numXcreated(ig,:) = [lonelyGroup(ig), num_breakup - 1];
        end
        
        % 2.6  In lonely group, find repetations (A245+A145, A145+A345, A245+A345)
        linearComb4 = linearComb3(lonelyGroup,:);
        for ifile = 1:N %
            K1 = find(demandType == char(files(ifile)));  % users requesting this file
            K2 = setxor(1:K, K1); % users not requesting this file
            if length(K1) >= t % this only happens when the # of users requesting this file >= t, 
                % happens any cases when W_{k1\in K1, k2\in K2}, |k1| <= t-1,
                for icandK1 = 1:t-1
                    candidate1 = int2str(K1(1:icandK1)); % (the first) icandK1 indices from K1
                    candidate2_all = nchoosek(K2, t - icandK1); % all possible choices from K2 
                    for icandK2 = 1:length(candidate2_all) % icandK2 indices from K2
                        candidate2 = int2str(candidate2_all(icandK2,:));
                        candidate2 = candidate2(find(~isspace(candidate2)));
                        candidate = strcat(char(files(ifile)), candidate1, candidate2); % A145, (t=2, (A,A,A,B,B,C))
                        % look for A145
                        if length(find(linearComb4 == candidate)) >= 2
                            A145Companion = nchoosek(K1, icandK1); % 1,2,3
                            rows = [];
                            for iComp = 1:length(A145Companion)
                                ind = find(linearComb4 == strcat(char(files(ifile)), int2str(A145Companion(iComp)), candidate2));
                                [row, ~] = ind2sub(size(linearComb4), ind);
                                rows = [rows; row];
                            end
                            rows = unique(rows); % all involving rows
                            [keepRow,~] = ind2sub(size(linearComb4), find(linearComb4 == candidate)); % A145+A245; A145+A345
                            duplicateRow = setxor(keepRow, rows); %A245+A345
                            % find which partition creates duplication transmission
                            aLinearComb = cellstr(linearComb3(lonelyGroup(duplicateRow(1)),:));
                            dispPartObj = func_waysToBreakUp(aLinearComb);
                            x_pos = func_xPos(dispPartObj, char(files(ifile)), t);
                            temp = zeros(size(numXcreated,1),length(x_pos));
                            temp(rows,:) = ones(length(rows),1)*x_pos';
                            numXcreated = [numXcreated, temp]; % 1st: linearComb, 2nd:each creates # of x, 3rd-..., which x correlaeted 
                        end
                    end
                end
            end
        end
%         for iX = 3:size(numXcreated, 2)
%             extraAeq = numXcreated(:, iX);
%             xind =  
%         end
        
        % 2.7  Find the joint groups (linear combinations), and AA, x, f:
        irow = 1;
        leftrow = 1:num_subsets;
        while 1
            posof1 = find(relationship(irow,:) == 1, 1); % if find one
            if ~isempty(posof1) % not empty
                group = irow;
                while 1
                    [all1i, all1j] = ind2sub([length(group), num_subsets], find(relationship(group,:) == 1)); % find all 1's in group row
                    if size(all1j, 1) > 1
                        all1j = all1j.';
                    end
                    relationship(group, :) = 0;
                    if isempty(all1j)
                        break;
                    else
                        group = [group, all1j]; % initialize group to be all the 1's in this row
                        group = unique(group); % exclude repetations
                    end
                end
                group  % show groups for debugging
                % ==finished============================================================
                aLinearComb = cellstr(linearComb3(group(1),:));
                dispPartObj = func_waysToBreakUp(aLinearComb);
                num_breakup = length(dispPartObj); % # of break-ups for this Y
                if num_breakup >= 2 
                    % ================================================================
                    for ig = 1:length(group) % each Y
                        aLinearComb = cellstr(linearComb3(group(ig),:));
                        dispPartObj = func_waysToBreakUp(aLinearComb);
                        temp = cellstr(linearComb3(group(ig),:));
                        temp = func_waysToBreakUp(temp);
                        storeItems{ig} = temp;
                        if ig == 1 % for the first Y in this group, do the same thing as lonelygroup
                            for ip = 2:num_breakup % for each breakup
                                breakup = dispPartObj{ip};
                                [a(ip-1,:), f(ip-1)] = func_genOptLonelyGroup(breakup, files, N, K, t);
                            end

                        else % for the other Y's, do not need to create new variable x, just change AA, F 

                            for ip = 2:num_breakup % for each breakup
                                breakup = dispPartObj{ip};
                                [a2(ip-1,:), f2(ip-1)] = func_genOptJointGroup(storeItems, breakup, files, N, K, t);
                            end
                            a = a + a2; f = f + f2;
                        end
                    end
                    aeq = (num_breakup > 2) * (num_breakup - 1); % more than one way to breakup, sum must be 1
                    AA = [AA; a]; F = [F, f]; x = [x; zeros(num_breakup - 1, 1)]; Aeq = [Aeq, aeq]; % for each breakup add a line to AA and f
                    clear a, clear a2, clear f, clear f2
                    % ================================================================
                end
                clear storeItems
                % ==finished============================================================
                leftrow = setxor(leftrow, group);
                if isempty(leftrow)
                    fprintf('There are no more ralationships.')
                    break
                else
                    irow = leftrow(1);
                end
            else
                leftrow = setxor(leftrow, irow);
                irow = irow + 1;
                if irow > num_subsets
                    fprintf('There are no more ralationships.')
                    break
                end
            end
        end
        
        %% 3  Save
        Aeq = func_aeq2aeq(Aeq, x);   
        demand = char(demandType);
        demand = demand(10:end-3); % delete 'matrix([[' and ']])'
        demand = strcat('(', demand, ')');
        filename = strcat('C:\Users\Kai\Dropbox\Caching\t=', int2str(t), ', demand=', demand, '.mat'); 
        save(filename,'AA','F','Aeq','x');
    end 
end




