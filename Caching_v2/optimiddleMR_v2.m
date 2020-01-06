clc; clear all; close all

%% 0.  Setting: 
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
% K = sym('K') % user
L = sym('L');
% M = sym('M'); % cache
% N = sym('N')  % file
O = sym('O');
P = sym('P');
Q = sym('Q');
% R = sym('R'); % rate
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
Ned = min(K,N); % only consider demand types where all files are requested
figure
% Tian
t_axis = 0:K;
M_Tian = t_axis.*((N-1)*t_axis+K-N)/K/(K-1);
R_Tian =  N*(K-t_axis)/K;
plot(M_Tian,R_Tian, 'k--d', 'Linewidth', 1) 

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
plot(M_Yu,R_Yu,'r-x', 'Linewidth', 1) 

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

%% 6.  Optimize (M, R)
t_axis = 0:K;
P_all = cell(length(t_axis)-2, 1);
for it = 2:length(t_axis) - 1    
    % 6.1  read in data from each demandType
    t = t_axis(it);
    M_ = []; R_ = []; Aeq_ = []; x_ = [];
    for ip = 1:size(partition, 1)
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
        M_ = blkdiag(M_, M); R_ = blkdiag(R_, R); Aeq_ = blkdiag(Aeq_, Aeq); x_ = [x_; x];
        clear M, clear R, clear Aeq, clear x
    end
    % 6.2  Initialize \script{P} to be two extreme points    
    R = R_/nchoosek(K,t); % normalize
    M = M_/nchoosek(K,t); % normalize
    P = [M_Yu(it), R_Yu(it); % initialize
         M_Tian(it), R_Tian(it)];
    Px = [zeros(size(x_, 1), 1), func_calcEndPointx(Aeq_, x_)]; % 's', 't' are missing here
    Aeq = [Aeq_, zeros(size(Aeq_,1),2)]; % additional variables 's' and 't'
    beq = ones(size(Aeq_,1), 1); % do not need add 's' and 't' in b
    x = [x_; 0; 0]; % additional variable 's' and 't'


    % 6.3  optimize
    n = 2; ii = 1; 
    while ii < n
        alpha = (P(ii+1,1) - P(ii,1))/(P(ii,2) - P(ii+1,2));
        beta = P(ii,1) + alpha*P(ii,2);
        % optimization
        F = zeros(size(x)); F(end-1:end) = [-1, alpha]; % min(- s + alpha*t)
        AA1 = -[M; -1*ones(1,ip*K); zeros(1,ip*K);]; % format: vertical
        AA2 = [R'; zeros(1, ip); -1*ones(1, ip)]; % format: vertical
        AA = [AA1, AA2]'; % transpose
        b = zeros(size(AA, 1),1); % Mx >= s, Rx <= t;
        lb = [zeros(size(x,1) - 2, 1); 0; 0]; % variable 's' at least 0, 't' at least 0,
        [ub, beq] = func_calcUbBeq(Px(:,ii), Px(:,ii+1), Aeq(:,1:end-2), beq);
        ub = [ub; P(ii,1)-P(ii+1,1); P(ii+1,2)-P(ii,2)]; % variable 's' at most adds deltaM, 't' at most adds deltaR 
        options = optimoptions('linprog','Algorithm','dual-simplex'); % Set options to use the 'dual-simplex' algorithm.
        [x, fval, exitflag, output] = linprog(F,AA,b,Aeq,beq,lb,ub)  %,options)  
        demandName = func_findBetaNew(M'*x(1:end-2), R*x(1:end-2), alpha, ip, K);
        betanew = P(ii,1) + alpha*P(ii,2) + alpha*R(demandName,:)*x(1:end-2) - min(M(:, 1+(demandName-1)*K: demandName*K)'*x(1:end-2));
        M_new = P(ii,1) - min(M(:, 1+(demandName-1)*K: demandName*K)'*x(1:end-2));
        R_new = P(ii,2) + R(demandName, :)*x(1:end-2);
        if (betanew < beta) & (abs(betanew-beta) > 0.001)
            P = [P(1:ii,:); [M_new, R_new]; P(ii+1:end,:)];
            Px = [Px(:,1:ii), Px(:,1:ii)+x(1:end-2), Px(:,ii+1:end)]; % Px does not include 's'
            n = n + 1;
%                 plot(M_new, R_new, 'k-d')
        else
            ii = ii + 1;
        end
    end
    P_all{it-1} = P;
    clear R, clear M, clear P, clear Px, clear Aeq, clear beq
end


for ip = 1:length(P_all)
    points = P_all{ip}; 
    if size(points,1) >= 3 % there are new points
        x = points(:,1); y = points(:,2);
        plot(x,y, 'b-o', 'Linewidth', 1)
    end
end
legend('Tian-Chen scheme','Yu et al. scheme','Partial decomposition' )
ylabel('R')
xlabel('M')

    

