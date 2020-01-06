% function plotMR(N, K, Ned)
clc; close all; clear all

%% 0.  Setting: 
N = 6; K = 9;
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
% K = sym('K') % user K
L = sym('L');
% M = sym('M'); % cache M
% N = sym('N')  % file N
O = sym('O');
P = sym('P');
Q = sym('Q');
% R = sym('R'); % rate R
S = sym('S');
T = sym('T');
U = sym('U');
V = sym('V');
W = sym('W');
X = sym('X');
Y = sym('Y');
Z = sym('Z');
files = [A,B,C,D,E,F,G,H,I,J,L,O,P,Q,S,T,U,V,W,X,Y,Z];

%% 1.  Plot Tian v.s. Yu:
Ned = min(K, N); % only consider demand types where all files are requested
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

%% 2.  Demand types
if K - N == 0
    partition = 0;
elseif K - N == 1
    partition = 1;
elseif K - N >= 2
    partition = intpartgen(K - N, N);
    partition = cell2mat(partition(end));
end
if size(partition, 2) < N % each file requested by at least one user, the rest (K-N) partition assigned to first users
    partition = [partition, zeros(size(partition,1), N - size(partition, 2))];
end
partition = ones(size(partition)) + partition; % all files are requested: each file at least one user

%% 3. Iterate each t (=2 for testing)
t_axis = 0:K;
for it = 2:length(t_axis) - 1 
    t = t_axis(it);
    %% 4. Iterate each demand (= (A,A,B,B,C,C) for testing):
    for ip = 1:size(partition, 1)
        AA = []; F = []; x = []; Aeq = [];
        demandType = [];
        for sym = 1:N
            demandType = [demandType, files(sym)*ones(1, partition(ip, sym))];
        end
        leader = zeros(1,N);
        for ifile = 1:N % leader for each file
            leader(ifile)= find(demandType == files(ifile), 1);
        end
        
        
        % 4.1.  t+1 transmissions (label of Y):
        t_plus_1_subsets = nchoosek(1:K, t+1); % all t+1 subsets
        % 4.2.  Delete redundant Y:
        if length(leader) > K
            fprintf('There are no redundant linear combinations')
        elseif t+1 > K - length(leader)
            fprintf('There are no redundant linear combinations')
        else
            non_leader = setxor(1:K, leader);
            licombDel = nchoosek(non_leader, t+1); % linear combinations ready to delete
            [t_plus_1_subsets, ~] = setdiff(t_plus_1_subsets, licombDel, 'rows');
        end
        
        
        % 4.3  Symbols in each Y:
        num_subsets = size(t_plus_1_subsets, 1); % # of Y
        linearComb = strings([num_subsets, t+1]); % initialize linear combination
        linearCombFile = linearComb; % This is a copy
        for iLiComb = 1:num_subsets % each Y
            t_plus_1_subset = t_plus_1_subsets(iLiComb,:);
            t_subset = nchoosek(t_plus_1_subset,t); % each symbol in this Y
%             % 4.3.1  Create symbols in this Y:
            for iSymbol = 1:length(t_subset) 
                filename = demandType(setdiff(t_plus_1_subset, t_subset(iSymbol,:)));
                symbol = strcat(char(filename),strjoin(string(t_subset(iSymbol,:)),''));
                linearComb(iLiComb, iSymbol) = symbol; % this is for finding groups, will be ruined
                linearCombFile(iLiComb, iSymbol) = uint8(char(filename)); % this is for checking singles, will be ruined
            end
        end
        
        % 4.4  Finding groups
        groups = func_findGroups(linearCombFile, t);
        num_variables = []; M = []; R = [];
        for ig = 1:length(groups)
            % 4.5  for each group, find breakups
            group = groups{ig}; % e.g. group = [5,6]
            groupPartation = cell(size(group));
            for ig = 1:length(group)
                aLinearComb = cellstr(linearComb(group(ig),:));
                temp = func_waysToBreakUp(aLinearComb);
                for it = 1:length(temp)
                    groupPartationArray{ig, it} = temp{it};
                end
            end
            
            % 4.6  create F(i.e. R), A(i.e. M)
            [num_variables_temp, M_temp, R_temp] = func_createX(groupPartationArray, t, demandType, leader);
            num_variables = [num_variables, num_variables_temp];
            M = [M; M_temp];
            R = [R, R_temp];
            clear num_variables_temp, clear M_temp, clear R_temp, clear groupPartationArray
        end

        % 4.7 calculate Aeq
        Aeq = func_calcAeq(num_variables, sum(num_variables));
        x = zeros(sum(num_variables), 1);
        
        
        %% 5.  Save
        demand = char(demandType);
        demand = demand(10:end-3); % delete 'matrix([[' and ']])'
        demand = strcat('(', demand, ')');
        filename = strcat('C:\Users\Kai\Dropbox\Caching\N=',int2str(N),',K=',int2str(K),',t=', int2str(t), ', demand=', demand, '.mat'); 
        save(filename,'M','R','Aeq','x');
    end 
end

%% 6.  Optimize (M, R)
t_axis = 0:K;
P_all = cell(length(t_axis)-2, size(partition, 1));
for it = 2:length(t_axis) - 1    
    % 6.1  read in data from each demandType
    t = t_axis(it);
    for ip = 1:size(partition, 1)
        M_ = []; R_ = []; Aeq_ = []; x_ = [];
        demandType = [];
        for sym = 1:N
            demandType = [demandType, files(sym)*ones(1, partition(ip, sym))]; 
        end
        demand = char(demandType);
        demand = demand(10:end-3); % delete 'matrix([[' and ']])'
        demand = strcat('(', demand, ')');
        filename = strcat('C:\Users\Kai\Dropbox\Caching\N=',int2str(N),',K=',int2str(K),',t=', int2str(t), ', demand=', demand, '.mat'); 
        load(filename);
        % ============================================
        M_ = [M_; M]; R_ = [R_, R]; Aeq_ = blkdiag(Aeq_, Aeq); x_ = [x_; x];
        clear M, clear R, clear Aeq, clear x
    
        % 6.2  Initialize \script{P} to be two extreme points    
        R = R_/nchoosek(K,t); % normalize
        M = M_/nchoosek(K,t); % normalize
        P = [M_Yu(it), R_Yu(it); % initialize
             M_Tian(it), R_Tian(it)];
        Px = [zeros(size(x_, 1), 1), func_calcEndPointx(Aeq_, x_)]; % 's' is missing here
        Aeq = [Aeq_, zeros(size(Aeq_,1),1)]; % additional variable 's'
        beq = ones(size(Aeq_,1), 1); % do not need add 's' in b
        x = [x_; 0]; % additional variable 's'


        % 6.3  optimize
        n = 2; ii = 1; 
        while ii < n
            alpha = (P(ii+1,1) - P(ii,1))/(P(ii,2) - P(ii+1,2));
            beta = P(ii,1) + alpha*P(ii,2);
            % optimization
            F = [alpha*R, -1]; % min(alpha*Rx - s)
            AA = -[M;-1*ones(1,K)]'; b = zeros(size(AA, 1),1); % min(Mx)>=s, i.e. M'x - ones(6,1)*s >= 0, [M, -1*ones(1,6)]'*[x, s]T >= 0
            lb = [zeros(size(x,1)-1,1); 0]; % variable 's' at least 0,
            [ub, beq] = func_calcUbBeq(Px(:,ii), Px(:,ii+1), Aeq(:,1:end-1), beq);
            ub = [ub; P(ii,1)-P(ii+1,1)]; % variable 's' at most adds deltaM  
            options = optimoptions('linprog','Algorithm','dual-simplex'); % Set options to use the 'dual-simplex' algorithm.
            [x, fval, exitflag, output] = linprog(F,AA,b,Aeq,beq,lb,ub);  %,options)    
            betanew = P(ii,1) + alpha*P(ii,2) + alpha*R*x(1:end-1) - min(M'*x(1:end-1));
            M_new = P(ii,1) - min(M'*x(1:end-1));
            R_new = P(ii,2) + R*x(1:end-1);
            if (betanew < beta) & (abs(betanew-beta) > 0.001)
                P = [P(1:ii,:); [M_new, R_new]; P(ii+1:end,:)];
                Px = [Px(:,1:ii), Px(:,1:ii)+x(1:end-1), Px(:,ii+1:end)]; % Px does not include 's'
                n = n + 1;
%                 plot(M_new, R_new, 'k-d')
            else
                ii = ii + 1;
            end
        end
        P_all{it-1, ip} = P;
        clear R, clear M, clear P, clear Px, clear Aeq, clear beq
    end
    func_plotMR(P_all(it-1,:), size(partition, 1))
end






