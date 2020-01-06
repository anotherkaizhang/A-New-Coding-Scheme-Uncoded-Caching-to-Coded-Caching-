clc; clear all; close all
%% 0.1  Plot Tian v.s. Yu:
N = 3; K = 4;
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
% 0.2 Demand types:
partition = intpartgen(K - N);
partition = cell2mat(partition(end));
if size(partition, 2) < N
    partition = [partition, zeros(size(partition,1), N - size(partition, 2))];
end
partition = ones(size(partition)) + partition; % all files are requested: each file at least one user

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

plot(M_Yu(4), R_Yu(4),'k-d', 'Linewidth', 3)
DeltaM = 0.01;
%% 2  optimize (M,R)
t_axis = 0:K;
for it = 4 %2:length(t_axis) - 1 
    start_point = [M_Yu(it), R_Yu(it)]; % M_Yu
    end_point = [M_Tian(it), R_Tian(it)]; % Tian
    t = t_axis(it);
    %% 2. Iterate each demand (= (A,A,B,B,C,C) for testing):
    for ip = 2% 1:size(partition, 1)
        demandType = [];
        for sym = 1:N
            demandType = [demandType, files(sym)*ones(1, partition(ip, sym))]; 
        end
        demand = char(demandType);
        demand = demand(10:end-3); % delete 'matrix([[' and ']])'
        demand = strcat('(', demand, ')');
        filename = strcat('C:\Users\Kai\Dropbox\Caching\t=', int2str(t), ', demand=', demand, '.mat'); 
        load(filename);
        % ============================================ 
        AA = - AA.'/nchoosek(K,t);
        b = -DeltaM*ones(K,1);       
        beq = ones(size(Aeq,1), 1);
        lb = zeros(1,length(x));
        ub = ones(1,length(x));
        F = -F/nchoosek(K,t);
        
        options = optimoptions('linprog','Algorithm','dual-simplex'); % Set options to use the 'dual-simplex' algorithm.
        while ~isempty(x)  
            %Solve the linear program and request the function value, exit flag, and output structure.
            [x, fval, exitflag, output] = linprog(F,AA,b,Aeq,beq,lb,ub,options)
            x = x > 0; % maximize trip
            plot([start_point(1), start_point(1) + max(AA*x)], [start_point(2), start_point(2) - F*x],'r-*' )
            % delete the x = 1 positions for the next iteration 
            x = func_findInd(x, Aeq); 
            x(find(x == 1)) = [];
            if isempty(x) % reach Tian point
                fprintf('Reach destination.')
                isbreak = 1;
                break
            end
            AA(:,ind) = []; Aeq(:,ind) = []; 
            if isempty(Aeq)
                beq = [];
            else
                zerolines = find(any(Aeq==0,2)==0);
                Aeq(zerolines,:) = []; beq(zerolines,:) = [];
            end   
            lb = zeros(1,length(x));
            ub = ones(1,length(x));
            F(ind) = [];
            start_point(1) = start_point(1) + min(AA*x);
            start_point(2) = start_point(2) + F*x;
        end
        % ============================================
        if isbreak == 1; isbreak = 1; break; end
    end
    if isbreak == 1; isbreak = 1; break; end
end



%% continue
% new_start_point = [start_point(1) + min(A*x), start_point(2) + f*x];
% DeltaM = (new_start_point(1) - end_point(1))/1.5; % minimum step
% A = -[1,1,1,1,0,0;
%     1,1,1,1,0,0;
%     1,0,1,0,1,1;
%     1,0,1,0,1,1;
%     0,1,0,1,1,1;
%     0,1,0,1,1,1]/15;
% b = -DeltaM*ones(6,1);
% Aeq = [];
% beq = [];
% lb = zeros(1,6);
% ub = ones(1,6);
% f = ones(1,6)/15;
% % Set options to use the 'dual-simplex' algorithm.
% options = optimoptions('linprog','Algorithm','dual-simplex');
% %Solve the linear program and request the function value, exit flag, and output structure.
% [x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,lb,ub,options)
% % x = x>0;
% % x = ones(6,1);
% plot(new_start_point(1) + min(A*x), new_start_point(2) + f*x,'r-*' )

%% t = 2 black line:
M = [15, 11]/15;
R = [19, 25]/15;
plot(M,R,'g-d', 'Linewidth', 1.1) 

M = [15, 11]/15;
R = [19, 24]/15;
plot(M,R,'g-d', 'Linewidth', 1.1) 

%% t = 3 green line:
M = [30, 28]/20;
R = [15, 18]/20;
plot(M,R,'g-d', 'Linewidth', 1.1) 

M = [30, 20]/20;
R = [15, 27]/20;
plot(M,R,'g-d', 'Linewidth', 1.1) 