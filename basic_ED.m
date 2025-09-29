%% 读入文件
mpcfile = 'case39_ED.m';
loadfile = 'load3996.mat';
mpc = loadcase(mpcfile);
load(loadfile);

%% 定义参数
T = size(Pd,2);
num_bus = size(mpc.bus,1);       
num_gen = size(mpc.gen,1);
num_branch = size(mpc.branch,1);
Pgmax =mpc.gen(:,9);
Pgmin = mpc.gen(:,10);
G = zeros(num_bus,num_gen);                                                % 发电机位置
for i=1:num_gen
    G(mpc.gen(i,1),i) = 1;
end
PTDF = makePTDF(mpc);
Plmax = mpc.branch(:,6);
Plmin = -Plmax; 
U = zeros(T,T-1);                                                          % 实现x(t)-x(t-1)
for i=1:size(U,2)
    U(i,i) = -1;
    U(i+1,i) = 1;
end
Ru = 0.1*Pgmax;                                                            % 爬坡率
Rd = Ru;                                                                   % 此处认为向上向下爬坡能力一样
f = zeros(num_bus,3);                                                      % 运行成本
f(locat_gen,3) = mpc.gencost(:,5);
f(locat_gen,2) = mpc.gencost(:,6);
f(locat_gen,1) = mpc.gencost(:,7);                

%% ED
% 定义决策变量
Pg = sdpvar(num_bus,T,'full');
Pl = sdpvar(num_branch,T,'full');

% 添加约束条件
Constraints = [];
Constraints = [Constraints ; (sum(Pg,1)==sum(Pd,1)):'Pb'];                                                          % 负荷平衡约束
Constraints = [Constraints ; repmat(Plmin,1,T)<=Pl<=repmat(Plmax,1,T) ; (Pl==PTDF*(Pg-Pd)):'Pl'];                   % 直流潮流约束，其实不要
Constraints = [Constraints ; repmat(Pgmin,1,T)<=Pg<=repmat(Pgmax,1,T)];                                             % 出力上下限
Constraints = [Constraints ; -repmat(Rd,1,T-1)<=Pg*U<=repmat(Ru,1,T-1)];                                            % 爬坡约束

% 定义目标函数
Objective = sum(f(:,2)'*Pg);                                               % 成本设置为线性

%% 求解
optimize(Constraints,Objective);
result.cost = value(Objective);
result.PG = value(Pg);
result.PL = value(Pl);

%% 画图
% 机组组合
figure
t = 1:T;
plot(t,result.PG(30:39,:))
legend('gen1','gen2','gen3','gen4','gen5','gen6','gen7','gen8','gen9','gen10');


