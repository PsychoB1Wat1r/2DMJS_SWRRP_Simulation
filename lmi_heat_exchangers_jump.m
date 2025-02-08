%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%% 数值求解代码 + 无额外变量 + 2024-4-3 %%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%% 2D MJSs + 滑动窗口轮询协议（无零阶保持器） + 异步 %%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 无额外变量 模态-协议-依赖 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
heat_exchangers_jump % 导入系统参数
%% 循环定义决策变量
% gamma = sdpvar(1, 1,'symmetric');   % 待求的 H∞ 性能指标 gamma
gamma = 2.1603^2;
for i = 1 : numSysMode
    eval(['Ph_i_', int2str(i), '= sdpvar(1, 1,''symmetric'');']);    % 李雅普诺夫矩阵 Ph_ip
    eval(['Pv_i_', int2str(i), '= sdpvar(1, 1,''symmetric'');']);    % 李雅普诺夫矩阵 Pv_ip
    for r = 1 : numConMode
        for n = 1 : numProMode 
            eval(['Omegah_in_', int2str(i), int2str(n) '= sdpvar(na/2, na/2,''symmetric'');']);  % 模态相关事件触发权值矩阵 Omegah_ipn
            eval(['Omegav_in_', int2str(i), int2str(n) '= sdpvar(na/2, na/2,''symmetric'');']);  % 模态相关事件触发权值矩阵 Omegav_ipn
            eval(['G_ir_', int2str(i), int2str(r),'= sdpvar(na, na,''symmetric'');']);   % 松弛矩阵 G_ipgr
            eval(['Wh_rn_', int2str(r), int2str(n),'= sdpvar(mb/2, mb/2,''full'');']); % 解耦矩阵 Wh_grn
            eval(['Wv_rn_', int2str(r), int2str(n),'= sdpvar(mb/2, mb/2,''full'');']); % 解耦矩阵 Wv_grn
            eval(['K_rn_', int2str(r), int2str(n),'= sdpvar(numProMode, mb,''full'');']); % 控制增益 K_grn
        end
    end
end
% 高级转移概率 条件概率1----QG  转移概率 条件概率2----IR  轮询期刊----N
%% 中间变量预定义
% Schur补 对角矩阵 P_i
P_i_1 = blkdiag(Ph_i_1, Pv_i_1);  P_i_2 = blkdiag(Ph_i_2, Pv_i_2);  daigP = blkdiag(P_i_1, P_i_2);
% Schur补 对角矩阵 G_ir
G_i_1 = blkdiag(G_ir_11, G_ir_12);  G_i_2 = blkdiag(G_ir_21, G_ir_22);
% %% 循环定义LMI_1
for i = 1 : numSysMode
        Ph_i = eval(['Ph_i_', int2str(i)]);
        Pv_i = eval(['Pv_i_', int2str(i)]);
        P_i = blkdiag(Ph_i, Pv_i);
        Comb_P = [sqrt(varepsilon(i,1))*P_i   sqrt(varepsilon(i,2))*P_i];
        G_i = eval(['G_i_', int2str(i)]);
        LMI_1 = [ -P_i        Comb_P;
                  (Comb_P)'    -G_i];
        eval(['LMI_1_', int2str(i), '= [LMI_1 <= delta*eye(size(LMI_1))];']); 
end
%% 循环定义LMI_2
for i = 1 : numSysMode
    A = eval(['A', int2str(i)]);  B = eval(['B', int2str(i)]);  E = eval(['E', int2str(i)]);
    C = eval(['C', int2str(i)]);  D = eval(['D', int2str(i)]);  F = eval(['F', int2str(i)]);
    Ph_i = eval(['Ph_i_', int2str(i)]);
    Pv_i = eval(['Pv_i_', int2str(i)]);
    P_i = blkdiag(Ph_i, Pv_i);
    for r = 1 : numConMode
        for n = 1 : numProMode 
            Omegah_in = eval(['Omegah_in_', int2str(i), int2str(n)]);
            Omegav_in = eval(['Omegav_in_', int2str(i), int2str(n)]);
            Omega_in = blkdiag(Omegah_in, Omegav_in);
            Active_n = eval(['Active_n_', int2str(n)]);
            G_ir = eval(['G_ir_', int2str(i), int2str(r)]);

            Wh_rn = eval(['Wh_rn_', int2str(r), int2str(n)]);
            Wv_rn = eval(['Wv_rn_', int2str(r), int2str(n)]);
            W_rn = blkdiag(Wh_rn, Wv_rn);

            K_rn = eval(['K_rn_', int2str(r), int2str(n)]);

            Comb = [A*W_rn + B*Active_n*K_rn   -B*Active_n*K_rn   E];    

            theta_1 = blkdiag(theta_h(i,1), theta_v(i,1));  theta_2 = blkdiag(theta_h(i,2), theta_v(i,2));    % 区间分段1 

            Pi_13 = [Comb'*sqrt(theta_1)   Comb'*sqrt(theta_2)]';
            Pi_23 = [C*W_rn + D*Active_n*K_rn    -D*Active_n*K_rn    F];
            Pi_33 = blkdiag(beta*Omega_in + G_ir - W_rn - W_rn',    -Omega_in,    -gamma*eye(1,1));

            LMI_2 = [ -daigP          zeros(4,1)   Pi_13;
                       zeros(1,4)    -eye(1,1)     Pi_23;
                      (Pi_13)'       (Pi_23)'      Pi_33];
            eval(['LMI_2_', int2str(i), int2str(r), int2str(n),'= [LMI_2 <= delta*eye(size(LMI_2)), P_i >= delta*eye(size(P_i)), G_ir >= delta*eye(size(G_ir)), Omega_in >= delta*eye(size(Omega_in))];']); 
        end
    end
end
%% 线性矩阵不等式拼接
lmi1_Str = ''; % 初始化空字符，用于拼接
for i = 1 : numSysMode
    str_temp = eval(['LMI_1_', int2str(i), ',']);
    lmi1_Str = [lmi1_Str, str_temp]; % 线性矩阵不等式拼接
end
lmi2_Str = ''; % 初始化空字符，用于拼接
for i = 1 : numSysMode
    for r = 1 : numConMode
        for n = 1 : numProMode 
            str_temp = eval(['LMI_2_', int2str(i), int2str(r), int2str(n), ',']);
            lmi2_Str = [lmi2_Str, str_temp]; % 线性矩阵不等式拼接
        end
    end
end
%% 线性约束条件，需要同时满足
Constraint = [gamma >= 0, lmi1_Str, lmi2_Str];
sdpsettings('solver', 'mosek', 'verbos', 0); % 设置求解器为 mosek，并打印少量信息
reuslt = optimize(Constraint, gamma); % 自动求解的优化函数
if reuslt.problem == 0 % problem = 0 代表求解成功
    [primal,~] = check(Constraint); % 检查约束
    if min(primal) >= 0 && all(primal([1:2,4]) > 0) % 判断残差值
        disp('***********************************');
        disp('****Constraints are guaranteed*****');
        disp('**Primal residual is/are positive**');
        disp('***********************************');
        disp('****H∞ performance index*****');
        gamma = sqrt(double(gamma))
    else
        disp('********************************************');
        disp('**Warning: Primal residual is/are negative**');
        disp('********************************************');
        gamma = sqrt(double(gamma))
%         check(Constraint); % 检查残余值
    end
else
    disp('**************************************');
    disp('****Constraints are not guaranteed****');
    disp('**************************************'); 
    reuslt.info
    yalmiperror(reuslt.problem) % 打印出错信息
%     check(Constraint); % 检查残余值
end
%% 求解控制增益
for r = 1 : numConMode
    for n = 1 : numProMode 
        Wh_rn = eval(['Wh_rn_', int2str(r), int2str(n)]);
        Wv_rn = eval(['Wv_rn_', int2str(r), int2str(n)]);
        W_rn = blkdiag(Wh_rn, Wv_rn);
        K_rn = eval(['K_rn_', int2str(r), int2str(n)]);
        eval(['K_rn_', int2str(r), int2str(n), '= value(K_rn)*inv(value(W_rn));']);
    end
end
K_rn_11(isnan(K_rn_11)) = 0;  K_rn_21(isnan(K_rn_21)) = 0;  
K_rn_12(isnan(K_rn_12)) = 0;  K_rn_22(isnan(K_rn_22)) = 0;
K_rn_13(isnan(K_rn_13)) = 0;  K_rn_23(isnan(K_rn_23)) = 0;
K_rn_14(isnan(K_rn_14)) = 0;  K_rn_24(isnan(K_rn_24)) = 0;
disp('**************************************');
disp('******节点 1 2 被激活时的控制增益K******');
disp('**************************************');
K_rn_11,  K_rn_21   % 节点 1 2 被激活
disp('**************************************');
disp('******节点 2 3 被激活时的控制增益K******');
disp('**************************************');
K_rn_12,  K_rn_22   % 节点 2 3 被激活
disp('**************************************');
disp('******节点 3 4 被激活时的控制增益K******');
disp('**************************************');
K_rn_13,  K_rn_23   % 节点 1 3 被激活
disp('**************************************');
disp('******节点 4 1 被激活时的控制增益K******');
disp('**************************************');
K_rn_14,  K_rn_24   % 节点 1 3 被激活
%% 求解事件触发矩阵
for i = 1 : numSysMode
    for n = 1 : numProMode 
        Wh_rn = eval(['Wh_rn_', int2str(i), int2str(n)]);
        Wv_rn = eval(['Wv_rn_', int2str(i), int2str(n)]);
        W_rn = blkdiag(Wh_rn, Wv_rn);
        trans_W_rn = W_rn';
        Omegah_in = eval(['Omegah_in_', int2str(i), int2str(n)]);
        Omegav_in = eval(['Omegav_in_', int2str(i), int2str(n)]);
        Omega_in = blkdiag(Omegah_in, Omegav_in);
        eval(['Omega_in_', int2str(i), int2str(n), '= inv(value(trans_W_rn))*value(Omega_in)*inv(value(W_rn)); ']);
    end
end
disp('**************************************');
disp('**节点 1 2 被激活时的事件触发矩阵Omega**');
disp('**************************************');
Omega_in_11,  Omega_in_21
disp('**************************************');
disp('**节点 2 3 被激活时的事件触发矩阵Omega**');
disp('**************************************');
Omega_in_12,  Omega_in_22
disp('**************************************');
disp('**节点 3 4 被激活时的事件触发矩阵Omega**');
disp('**************************************');
Omega_in_13,  Omega_in_23
disp('**************************************');
disp('**节点 4 1 被激活时的事件触发矩阵Omega**');
disp('**************************************');
Omega_in_14,  Omega_in_24
disp(['运行时间: ', num2str(toc)]); 