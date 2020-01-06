% function plotMR(N, K, Ned)
clc; close all; clear all

%% 0. Setting: 
N = 3; K = 6;
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
for it = 2% 2:length(t_axis) - 1 
    start_point = [M_Yu(it), R_Yu(it)]; % M_Yu
    end_point = [M_Tian(it), R_Tian(it)]; % Tian
    t = t_axis(it);
    %% 2. Iterate each demand (= (A,A,B,B,C,C) for testing):
    for ip = 2% 1:size(partition, 1)
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
            end
        end
                
        
        % 2.6  Find the joint groups (linear combinations), and AA, x, f:
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
                aLinearComb = cellstr(linearComb3(group(1),:));
                dispPartObj = func_waysToBreakUp(aLinearComb);
                num_breakup = length(dispPartObj) % # of break-ups for this Y % each cell is a way of partition
                if num_breakup > 1 % do nothing if there is only one breakup (do not break)
                    a = zeros(num_breakup - 1, K); f = zeros(1, num_breakup - 1);
                    for ip = 2:num_breakup
                        breakup = dispPartObj{ip};
                        [a(ip-1,:), f(ip-1)] = func_genOptLonelyGroup(breakup, files, N, K, t);
                    end
                    aeq = (num_breakup > 2) * (num_breakup - 1); % more than one way to breakup, sum must be 1
                    AA = [AA; a]; F = [F, f]; x = [x; zeros(num_breakup - 1, 1)]; Aeq = [Aeq, aeq]; 
                end
%                     % ===========================================================
%                     a = zeros(num_breakup - 1, K); f = zeros(1, num_breakup - 1);
%                     for ig = 1:length(group)
%                         aLinearComb = cellstr(linearComb3(group(ig),:));
%                         dispPartObj{ig} = func_waysToBreakUp(aLinearComb); % each cell is a way of partition              
%                         for ip = 2 : num_breakup
%                             if ig == 1
%                                 PreviousBreakups = [];
%                             else
%                                 for iprev = 1:ig - 1
%                                     PreviousBreakups{iprev} = dispPartObj{iprev}{ip};
%                                 end
%                             end
%                             breakup = dispPartObj{ig}{ip};
%                             [a(ip-1,:), f(ip-1)] = func_genOptJointGroup(PreviousBreakups, breakup, files, N, K, t);
%                         end
%                         aeq = (num_breakup > 2) * (num_breakup - 1); % more than one way to breakup, sum must be 1
%                         AA = [AA; a]; F = [F, f]; x = [x; zeros(num_breakup - 1, 1)]; Aeq = [Aeq, aeq];
%                     % ===========================================================
%                     end
%                     leftrow = setxor(leftrow, group);
%                     if isempty(leftrow)
%                         fprintf('There are no more ralationships.')
%                         break
%                     else
%                         irow = leftrow(1);
%                     end
                % 2.6.0  Find which symbol repeats in this group:
%                 [~,idx_last,~] = unique(linearComb3(group,:));
%                 repeat_symbol = linearComb3(group,:);
%                 repeat_symbol = strcat("{",repeat_symbol(setxor(1:length(group)*(t+1), idx_last)),"}");
%                 repeat_symbol = repeat_symbol{1};
                % 2.6.1  Each Y in this group:
                % 2.6.1  In each Y, put the same file together, and find all possible ways to break up:
%                 aLinearComb = ;
%                         for ip = 2:num_breakup % for each breakup
%                             if ig == 1 % the first Y in the group, there is no previous one
%                                 PreviousBreakups == [];
%                             else % the following Y's, 
%                                 PreviousBreakups = cell(ig-1);
%                                 for iprev = 1:ig-1
%                                     PreviousBreakups{iprev} = dispPartObj{iprev}{ip};
%                                 end
%                                 breakup = dispPartObj{ig}(ip);
%                             end
%                             [a(ip-1,:), f(ip-1)] = func_genOptJointGroup(PreviousBreakups, breakup, files, N, K, t);
%                             end
%                             aeq = (num_breakup > 2) * (num_breakup - 1); % more than one way to breakup, sum must be 1
%                             AA = [AA; a]; F = [F, f]; x = [x; zeros(num_breakup - 1, 1)]; Aeq = [Aeq, aeq]; % for each breakup add a line to AA and f
% 
%                                                     
% %                         if ig == 1 % for the first Y in this group, do the same thing as lonelygroup
% %                             aeq = (num_breakup > 2) * (num_breakup - 1); % more than one way to breakup, sum must be 1
% %                             AA = [AA; a]; F = [F, f]; x = [x; zeros(num_breakup - 1, 1)]; Aeq = [Aeq, aeq]; % for each breakup add a line to AA and f
% %                             % 'a' and 'f' are a little due to repeated symbol, decrease a little for the other Y's (ig = 2,3,...), precompute the following variables 
% %                             u = zeros(1,t);
% %                             for irept = 1:t
% %                                 u(irept) = str2num(repeat_symbol(irept + 2));
% %                             end
% %                             repeatLine = strfind(dispPartObj, repeat_symbol); % where repeat symbol appears
% %                             AAFdecreaseLine = zeros(1, length(dispPartObj));
% %                             for id = 1 : length(dispPartObj)
% %                                 AAFdecreaseLine(id) = (1 - isempty(repeatLine{id})) * (id-1);
% %                             end
% %                             AAFdecreaseLine = AAFdecreaseLine(AAFdecreaseLine~=0);
% %                         else % for the other Y's, do not need to create new variable x, just change AA, F 
% %                             a(AAFdecreaseLine, u) = a(AAFdecreaseLine, u) - 1;
% %                             f(AAFdecreaseLine) = f(AAFdecreaseLine) - 1;
% %                             AA(end - size(a,1) + 1:end, :) = AA(end - size(a,1) + 1:end, :) + a;
% %                             F(end - length(f) + 1: end) = F(end - length(f) + 1: end) + f; 
% %                         end
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
%         demand = char(demandType);
%         demand = demand(10:end-3); % delete 'matrix([[' and ']])'
%         demand = strcat('(', demand, ')');
%         filename = strcat('C:\Users\Kai\Dropbox\Caching\t=', int2str(t), ', demand=', demand, '.mat'); 
%         save(filename,'AA','F','Aeq','x');
    end 
end
% 
% DeltaM = 0.2; % minimum step
% plot(M_Yu(3), R_Yu(3),'k-d', 'Linewidth', 3)
% A = -[1,1,1,1,0,0,4,2,2,4;
%       1,1,1,1,0,0,3,2,2,4;
%       1,0,1,0,1,1,2,4,2,4;
%       1,0,1,0,1,1,2,3,2,4;
%       0,1,0,1,1,1,2,2,4,4;
%       0,1,0,1,1,1,2,2,3,4]/15;
% A = [A;0,0,0,0,0,0,1,1,1,1];
% b = [-DeltaM*ones(6,1);1];
% Aeq = [];
% beq = [];
% lb = zeros(1,10);
% ub = ones(1,10);
% f = [ones(1,6), 4*ones(1,3), 5]/15;
% % Set options to use the 'dual-simplex' algorithm.
% options = optimoptions('linprog','Algorithm','dual-simplex');
% %Solve the linear program and request the function value, exit flag, and output structure.
% [x, fval, exitflag, output] = linprog(f,A,b,Aeq,beq,lb,ub,options)
% x = x > 0; % maximize trip
% plot(start_point(1) + min(A*x), start_point(2) + f*x,'r-*' )
% 
% %% continue
% % new_start_point = [start_point(1) + min(A*x), start_point(2) + f*x];
% % DeltaM = (new_start_point(1) - end_point(1))/1.5; % minimum step
% % A = -[1,1,1,1,0,0;
% %     1,1,1,1,0,0;
% %     1,0,1,0,1,1;
% %     1,0,1,0,1,1;
% %     0,1,0,1,1,1;
% %     0,1,0,1,1,1]/15;
% % b = -DeltaM*ones(6,1);
% % Aeq = [];
% % beq = [];
% % lb = zeros(1,6);
% % ub = ones(1,6);
% % f = ones(1,6)/15;
% % % Set options to use the 'dual-simplex' algorithm.
% % options = optimoptions('linprog','Algorithm','dual-simplex');
% % %Solve the linear program and request the function value, exit flag, and output structure.
% % [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options)
% % % x = x>0;
% % % x = ones(6,1);
% % plot(new_start_point(1) + min(A*x), new_start_point(2) + f*x,'r-*' )
% 
% %% t = 2 black line:
% M = [15, 11]/15;
% R = [19, 25]/15;
% plot(M,R,'g-d', 'Linewidth', 1.1) 
% 
% M = [15, 11]/15;
% R = [19, 24]/15;
% plot(M,R,'g-d', 'Linewidth', 1.1) 
% 
% %% t = 3 green line:
% M = [30, 28]/20;
% R = [15, 18]/20;
% plot(M,R,'g-d', 'Linewidth', 1.1) 
% 
% M = [30, 20]/20;
% R = [15, 27]/20;
% plot(M,R,'g-d', 'Linewidth', 1.1) 



