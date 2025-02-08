%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% 系统参数 + 2024-4-3 %%%%%%%%%%%%%%%%%%%%%%%%%%
%% 数据来源：Le Van Hien  2017 IEEE TCASII （热交换机模型） 稍作修改：迎合马氏链 %%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 模态1
a1 = 1.2;  b11 = 0.7;  b12 = 0.6;   b13 = -0.7;  b14 = 0.55;
a2 = 0.6;  b21 = 0.5;  b22 = 0.55;  b23 = -0.65; b24 = 0.45;
dx = 0.1;  dt = 0.1;    % 有限差分步长
A1 = [0         1;
      dt/dx     1 - dt/dx - a1*dt];
B1 = [0  b11*dt;
      0  b12*dt;
      0  b13*dt;
      0  b14*dt]';  % 假设有三个执行器节点
E1 = [1.2  -1.5]';
C1 = [0     0.9];
D1 = [-0.1, -0.2, -0.3, 0.1];
F1 = -0.1;
% 模态2
A2 = [0         1;
      dt/dx     1 - dt/dx - a2*dt];
B2 = [0  b21*dt;
      0  b22*dt;
      0  b23*dt;
      0  b24*dt]';
D2 = D1;  E2 = E1;  C2 = C1;  F2 = F1;
%% 齐次马尔可夫链 转移概率 
theta_h = [0.2 0.8;  0.7 0.3];   % low-level 区间分段1 水平转移概率 (theta_h_1)
theta_v = [0.8 0.2;  0.3 0.7];    % low-level 区间分段1 垂直转移概率 (theta_v_1)
%% 隐马模型 条件转移概率
varepsilon = [0.8 0.2;  0.15 0.85];   % 异步行为 2 （参数 zeta 与 tau）
% varepsilon = [1 0;  0 1];   % 异步行为 2 （参数 zeta 与 tau）
% varepsilon = [1 0;  1 0];   % 异步行为 2 （参数 zeta 与 tau）
%% 双向事件触发参数
cnt = 0.2;  etch = cnt;  etcv = cnt; % 水平与垂直方向上的触发阈值
beta = blkdiag(etch, etcv);  
%% RRP 协议信道使用者决策矩阵
Active_n_1 = blkdiag(1, 1, 0, 0); % 节点 1 2 被激活
Active_n_2 = blkdiag(0, 1, 1, 0); % 节点 2 3 被激活
Active_n_3 = blkdiag(0, 0, 1, 1); % 节点 3 4 被激活
Active_n_4 = blkdiag(1, 0, 0, 1); % 节点 4 1 被激活
numProMode = 4;
%% 系统矩阵维度
[ma,na] = size(A1);  [mb,nb] = size(B1);  [mc,nc] = size(C1);   % 系统状态 Xmn
[me,ne] = size(E1);  [md,nd] = size(D1);  [mf,nf] = size(F1);   % 被控输出 Ymn
%% 下标值
[numSysMode, numSysMode] = size(theta_h);   % 系统模态个数
[numConMode, numConMode] = size(varepsilon);    % 控制器模态个数
%% 病态求解法
delta = 0;  % 正定编程参数