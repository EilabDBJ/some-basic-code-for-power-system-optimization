%% 读入文件
mpcfile = 'case39_UC.m';
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
Su = 0.3*Pgmax;                                                            % 启动最大升出力
Rd = Ru;                                                                   % 此处认为向上向下爬坡能力一样
Sd = Su;
f = zeros(num_gen,3);                                                      % 运行成本 此处用二次函数
f(:,3) = mpc.gencost(:,5);
f(:,2) = mpc.gencost(:,6);
f(:,1) = mpc.gencost(:,7);
fs = mpc.gencost(:,2:3);                                                   % 启停成本
To = 16;                                                                   % 启停时间
Tc = 16;
Uo = zeros(T,T-To);                                                        % 用于启停约束
for i=1:size(Uo,2)
    Uo(i+1:(i+To),i) = 1;
end
Uc = ones(T,T-Tc);
for i=1:size(Uc,2)
    Uc(i:(i+Tc),i) = 1;
end                                     

%% UC
% 定义决策变量
Pg = sdpvar(num_gen,T,'full');                                             % 发电机出力
Sg = binvar(num_gen,T,'full');                                             % 启停状态
Pl = sdpvar(num_branch,T,'full');                                          % 线路功率
Cs = sdpvar(num_gen,T-1,'full');                                           % 启停成本

% 添加约束条件
Constraints = [];
Constraints = [Constraints ; (sum(Pg,1)==sum(Pd,1)):'Pb'];                                                          % 负荷平衡约束
% Constraints = [Constraints ; repmat(Plmin,1,T)<=Pl ; Pl<=repmat(Plmax,1,T) ; (Pl==PTDF*(G*Pg-Pd)):'Pl'];            % 直流潮流约束
% Constraints = [Constraints ; repmat(Pgmin,1,T).*Sg<=Pg ; Pg<=repmat(Pgmax,1,T).*Sg];                                % 出力上下限
Constraints = [Constraints ; Pg*U<=repmat(Ru-Su,1,(T-1)).*Sg(:,1:(T-1))+repmat(Su,1,(T-1))];                        % 爬坡约束
Constraints = [Constraints ; -Pg*U<=repmat(Rd-Sd,1,(T-1)).*Sg(:,2:T)+repmat(Sd,1,(T-1))];
Constraints = [Constraints ; Cs>=repmat(fs(:,1),1,T-1).*(Sg*U) ; Cs>=-repmat(fs(:,2),1,T-1).*(Sg*U) ; Cs>=0];       % 启停成本，这里认为开启和关停成本一致
Constraints = [Constraints ; Sg*Uo>=To*Sg*U(:,1:(T-To)) ; (1-Sg)*Uc>=-Tc*Sg*U(:,1:(T-Tc))];                         % 启停约束

% 定义目标函数
% Objective = trace(Pg'*diag(f(:,3))*Pg) + sum(f(:,2)'*Pg) + sum(f(:,1)'*Sg) + sum(sum(Cs));         % 二次带常数项
Objective = (sum(f(:,2)'*Pg) + sum(f(:,1)'*Sg))/4 + sum(sum(Cs));                                    % 一次带常数项

%% 求解
optimize(Constraints,Objective);
result.cost = value(Objective);
result.PG = value(Pg);
result.SG = value(Sg);
result.PL = value(Pl);

