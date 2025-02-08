%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% 闭环系统绘图代码 + 2024-4-3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% 2D MJSs + 滑动窗口轮询协议（无零阶保持器） + 异步 %%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmi_heat_exchangers_jump; % 直接调用脚本 lmi_heat_exchangers_jump 中的数据，即控制增益与事件触发矩阵
tic  % 计时器
%% 绘图坐标 绘图线条 坐标轴字体 大小
fontSizeXY = 10;  lineWidth = 1;  fontSizeAxis = 15; markerSize = 1;  rand('state', 50);  randn('state', 50);
%% 触发次数初始化/数据减少率计算
triggerTimes_h = 0;  triggerTimes_v = 0; % 初始化触发次数
sumTriggerTimes_h = 0;  sumTriggerTimes_v = 0; % 初始化重复模拟时的触发次数，用于累加
numOfRepeat = 1; % 重复模拟次数
%% 绘图参数
setp = 40;  % 步长
T = setp;  dt = 1;  Nt = T/dt;  % 参数: 水平方向
Lx = setp;  dx = 1;  Nx = Lx/dx;	% 参数: 垂直方向
xh_array = zeros(Nt+1, Nx+1);   % 系统水平状态信号存储数组
xv_array = zeros(Nt+1, Nx+1);   % 系统垂直状态信号存储数组
w_array = zeros(Nt+1, Nx+1);    % 噪声信号存储数组
u1_array = zeros(Nt+1, Nx+1);    % 控制信号存储数组
u2_array = zeros(Nt+1, Nx+1);    % 控制信号存储数组
u3_array = zeros(Nt+1, Nx+1);    % 控制信号存储数组
u4_array = zeros(Nt+1, Nx+1);    % 控制信号存储数组
y_array = zeros(Nt+1, Nx+1);   % 控制输出信号存储数组
eventh_array = ones(Nt+1, Nx+1);   % 水平方向触发瞬时存储数组
eventv_array = ones(Nt+1, Nx+1);   % 垂直方向触发瞬时存储数组
eventh_time_array = zeros(Nt+1, Nx+1);    % 保存水平事件时间以计算采样时间
eventv_time_array = zeros(Nt+1, Nx+1);    % 保存水平事件时间以计算采样时间
y_sum = zeros(Nt+1, Nx+1); % H∞ 性能指标演化存储数组
xh_new = zeros(1, 1);  xv_new = zeros(1, 1); % 触发控制初始时刻系统状态值选取为零 current
xh_old = zeros(1, 1);  xv_old = zeros(1, 1); % 触发控制初始时刻系统状态值选取为零 last
SysSwitSeq = ones(Nt+1, Nx+1);   % 系统模式存储数组
ConSwitSeq = ones(Nt+1, Nx+1); % 控制器模式存储数组
%% 系统边界条件  2019 Neurocomputing Dehao Li Jinling Liang Fan Wang
xh_array(1:Nt, 1) = 4.2*cos(2:Nt+1)'*sin(1);
xv_array(1:Nt, 1) = 4.2*cos(2:Nt+1)'*sin(2);
xh_array(1, 1:Nx) = 4.2*tanh(1)*sin(2:Nx+1)';
xv_array(1, 1:Nx) = 4.2*tanh(2)*sin(2:Nx+1)';
%% 噪声信号
for i = 1 : Nt
    for j = 1 : Nx
        w_array(i, j) = -7/((i+1)*(j+1));
%         w_array(i,j) = sin(0.15*(i + j))*exp(-0.1*(i + j));
    end
end
%% Roesser-type 模型双向切换序列生成
flagtemp = 1; % 初始化 零时刻模态为 1
for i = 1 : Nt+1
    for j = 1 : Nx+1
        if flagtemp == 1 % 若当前时刻系统为模态 1 
            a_h = rand;  a_v = rand;
            if a_h < theta_h(1,1)        
                flagSys_h = 1;       % System jump to mode 1
            else    
                flagSys_h = 2;       % System jump to mode 2
            end
            if a_v < theta_v(1,1)        
                flagSys_v = 1;       % System jump to mode 1
            else    
                flagSys_v = 2;       % System jump to mode 2
            end
        else
            if a_h < theta(2,1)        
                flagSys_h = 1;       % System jump to mode 1
            else    
                flagSys_h = 2;       % System jump to mode 2
            end
            if a_v < theta(2,1)        
                flagSys_v = 1;       % System jump to mode 1
            else    
                flagSys_v = 2;       % System jump to mode 2
            end
        end
        SysSwitSeq(i+1,j) = flagSys_h;  SysSwitSeq(i,j+1) = flagSys_v; 
    end
end
%% 异步切换 序列生成
for i = 1 : Nt+1
    for j = 1 : Nx+1
        if SysSwitSeq(i,j) == 1 % 若当前时刻系统为模态 1 
            Asy = rand;
            if Asy < varepsilon(1,1)
                flagConAsy2 = 1;
            else
                flagConAsy2 = 2;
            end
        else
            if Asy < varepsilon(2,1)
                flagConAsy2 = 1;
            else
                flagConAsy2 = 2;
            end
        end
        ConSwitSeq(i,j) = flagConAsy2;
    end
end
%% 轮询协议更新规则 
rrpSeq = zeros(1, Nt+1); % 用于存储节点执行序列
rrpSeq(1, 1:Nt+1) = rem(0:Nt, numProMode) + 1; % rem 函数用于取余操作
%% 节点接入次序的存储序列
nodes_Seq = zeros(4, Nt+1); % 用于存储节点执行序列
%% 绘图/数据包减少率计算
for r = 1 : numOfRepeat % 重复模拟
    triggerTimes_h = 0;  triggerTimes_v = 0; % 触发次数计数器
    for i = 1 : Nt
        last_eventh = 0;  last_eventv = 0; % 水平与垂直方向上 上一次事件出现的时刻
        for j = 1 : Nx
            switch rrpSeq(1,j)  
                case 1
                    Active = Active_n_1;
                    K_r_1 = K_rn_11;  K_r_2 = K_rn_21;
                    Omega_i_1 = Omega_in_11;  Omega_i_2 = Omega_in_21;
                case 2
                    Active = Active_n_2;
                    K_r_1 = K_rn_12;  K_r_2 = K_rn_22;
                    Omega_i_1 = Omega_in_12;  Omega_i_2 = Omega_in_22;
                case 3
                    Active = Active_n_3;
                    K_r_1 = K_rn_13;  K_r_2 = K_rn_23;
                    Omega_i_1 = Omega_in_13;  Omega_i_2 = Omega_in_23;
                otherwise
                    Active = Active_n_4;
                    K_r_1 = K_rn_14;  K_r_2 = K_rn_24;
                    Omega_i_1 = Omega_in_14;  Omega_i_2 = Omega_in_24;
            end
            nodes_Seq(:,j) = diag(Active);
            if  ConSwitSeq(i,j) == 1 
                K = K_r_1;
            else
                K = K_r_2;
            end
            if SysSwitSeq(i,j) == 1  
                omegaMatrix = Omega_i_1;
                As = A1;  Bs = B1;  Es = E1;  Cs = C1;  Ds = D1;  Fs = F1;
            else
                omegaMatrix = Omega_i_2;
                As = A2;  Bs = B2;  Es = E2;  Cs = C2;  Ds = D2;  Fs = F2;
            end
            xh_error = xh_array(i,j) - xh_old;  xv_error = xv_array(i,j) - xv_old; % 初始化状态误差
          %% 事件触发生成器
            % 水平方向事件触发生成器
            if ( omegaMatrix(1,1)*norm(xh_error) > beta(1,1)*omegaMatrix(1,1)*norm(xh_array(i,j)) )   % 水平触发条件成立
                eventh_array(i,j) = 1;  % 用于绘制水平触发时刻分布图
                eventh_time_array(i,j) = j - last_eventh;  % 用于绘制触发时刻的分布图
                last_eventh = j;
                triggerTimes_h = triggerTimes_h + 1; % 水平触发次数加 1
            else  
                eventh_array(i,j) = 0; % 用于绘制水平触发时刻分布图
                eventh_time_array(i,j) = 0;
            end
            % 垂直方向事件触发生成器
            if ( omegaMatrix(2,2)*norm(xv_error) > beta(2,2)*omegaMatrix(2,2)*norm(xv_array(i,j)) )   % 垂直触发条件成立
                eventv_array(i,j) = 1; % 用于绘制垂直触发时刻分布图
                eventv_time_array(i,j) = j - last_eventv;  % 用于绘制触发时刻的分布图
                last_eventv = j;
                triggerTimes_v = triggerTimes_v + 1; % 垂直触发次数加 1
            else  
                eventv_array(i,j) = 0; % 用于绘制垂直触发时刻分布图
                eventv_time_array(i,j) = 0;
            end
          %% 测量信号更新
            if eventh_array(i,j) == 1
                xh_new = xh_array(i,j);
                xh_old = xh_new; % 最后一次触发时刻的系统状态更新
            else
                xh_new = xh_old;
            end
            if eventv_array(i,j) == 1
                xv_new = xv_array(i,j);
                xv_old = xv_new; % 最后一次触发时刻的系统状态更新
            else
                xv_new = xv_old;
            end
            u_array = Active*K*[xh_new  xv_new]'; % 控制信号更新
            u1_array(i,j) = u_array(1,1); % 用于绘制控制信号1状态演化轨迹
            u2_array(i,j) = u_array(2,1); % 用于绘制控制信号2状态演化轨迹
            u3_array(i,j) = u_array(3,1); % 用于绘制控制信号3状态演化轨迹
            u4_array(i,j) = u_array(4,1); % 用于绘制控制信号4状态演化轨迹
          %% 系统响应
            xh_array(i+1,j) = As(1,1)*xh_array(i,j) + As(1,2)*xv_array(i,j) + Bs(1,:)*u_array + Es(1,1)*w_array(i,j);
            xv_array(i,j+1) = As(2,1)*xh_array(i,j) + As(2,2)*xv_array(i,j) + Bs(2,:)*u_array + Es(2,1)*w_array(i,j);
            y_array(i,j) = Cs(1,1)*xh_array(i,j) + Cs(1,2)*xv_array(i,j) + Ds(1,:)*u_array + Fs*w_array(i,j);
            y_sum(i,j) = sqrt(norm(y_array(1:i,1:j), 2))/sqrt(norm(w_array(1:i,1:j), 2)); % y用于绘制性能指标函数图像（范数）
        end 
    end
    sumTriggerTimes_h = sumTriggerTimes_h + triggerTimes_h; % 水平触发次数 均值求解 暂存变量
    sumTriggerTimes_v = sumTriggerTimes_v + triggerTimes_v; % 垂直触发次数 均值求解 暂存变量
end
%% 水平方向系统闭环演化图
figure('name', '水平方向系统闭环演化图')
[x,y] = meshgrid(0:dt:T-1, 0:dx:Lx-1);
z = xh_array(1:dt:Nt, 1:dx:Nx);
s1 = surf(x, y, z, 'LineWidth', lineWidth, 'FaceAlpha', 1);  colormap(copper);
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$n$', 'Interpreter', 'latex','Fontsize', fontSizeAxis);
ylabel('$m$', 'Interpreter', 'latex','Fontsize', fontSizeAxis);
zlabel('$x^{h}_{m,n}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
axis([0 T 0 Lx]);  grid on;  grid minor;  box on;
%% 垂直方向系统闭环演化图
figure('name', '垂直方向系统闭环演化图')
[x,y] = meshgrid(0:dt:T-1, 0:dx:Lx-1);
z = xv_array(1:dt:Nt, 1:dx:Nx);
s2 = surf(x, y, z, 'LineWidth', lineWidth, 'FaceAlpha', 1);  colormap(copper);
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('$x^{v}_{m,n}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
axis([0 T 0 Lx]);  grid on;  grid minor;  box on;
%% 事件触发释放时刻 xh
figure('name', '事件触发释放时刻 xh');
h = stem3((0:dt:T-1), (0:dt:T-1), eventh_time_array(1:dt:T, 1:dt:Lx), 'linewidth', lineWidth/2, 'marker', '.', 'color', '#1C1678', 'markersize', 15);
xdata = get(h,'xdata');  ydata = get(h,'ydata');  zdata = get(h,'zdata');   % 通过句柄获取绘图数据的 x y z 坐标值
zdata(zdata == 0) = NaN;    % 将 z = 0 的值重新置为 NaN
set(h,'xdata', xdata,'ydata', ydata, 'zdata', zdata);     % 更新绘图数据 (x 和 y 坐标值)
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('释放间隔', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
grid on;  grid minor;  box on;
%% 事件触发释放时刻 xv
figure('name', '事件触发释放时刻 xv');
h = stem3((0:dt:T-1), (0:dt:T-1), eventv_time_array(1:dt:T, 1:dt:Lx), 'linewidth', lineWidth/2, 'marker', '.', 'color', '#1C1678', 'markersize', 15);
xdata = get(h,'xdata');  ydata = get(h,'ydata');  zdata = get(h,'zdata');   % 通过句柄获取绘图数据的 x y z 坐标值
zdata(zdata == 0) = NaN;    % 将 z = 0 的值重新置为 NaN
set(h,'xdata', xdata,'ydata', ydata, 'zdata', zdata);     % 更新绘图数据 (x 和 y 坐标值)
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('释放间隔', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
grid on;  grid minor;  box on;
%% 系统模式切换
figure('name', '系统模式切换')
[x,y] = meshgrid(0:dt:T, 0:dx:Lx);
z = SysSwitSeq(1:dt:T+1, 1:dt:Lx+1);
colormap(copper);  surf(x,y,z);  view(0,90);
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
axis([0 T 0 Lx 0.5 2.5]);
%% 异步切换信号
figure('name', '异步切换')
[x,y] = meshgrid(0:dt:T, 0:dx:Lx);
z = ConSwitSeq(1:dt:T+1, 1:dt:Lx+1); 
colormap(copper);  surf(x,y,z);  view(0,90);
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
axis([0 T 0 Lx 0.5 2.5]);
%% 节点接入时序图
figure('name', '节点接入时序图');
stem((0:dt:T), nodes_Seq(1, 1:Nt+1), 'filled', 'marker', 'd', 'color', '#000000', 'linewidth', lineWidth*1.3, 'linestyle', 'none', 'markersize', 4.5*markerSize);
hold on;
stem((0:dt:T), 2*nodes_Seq(2, 1:Nt+1), 'filled', 'marker', 'd', 'color', '#E72929', 'linewidth', lineWidth*1.3, 'linestyle', 'none', 'markersize', 4.5*markerSize);
hold on;
stem((0:dt:T), 3*nodes_Seq(3, 1:Nt+1), 'filled', 'marker', 'd', 'color', '#362FD9', 'linewidth', lineWidth*1.3, 'linestyle', 'none', 'markersize', 4.5*markerSize);
hold on;
stem((0:dt:T), 4*nodes_Seq(4, 1:Nt+1), 'filled', 'marker', 'd', 'color', '#E178C5', 'linewidth', lineWidth*1.3, 'linestyle', 'none', 'markersize', 4.5*markerSize);
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$k$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
h2 = legend('节点 1', '节点 2', '节点 3', '节点 4', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
set(h2,'Orientation', 'horizon', 'Box', 'on');
ylabel('令牌分配', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
axis([0 T 0.5 5]);  grid on;  grid minor;  box on;
%% y 性能比较
figure('name', 'y 性能比较');
[x,y] = meshgrid(1:dt:T, 1:dx:Lx);  z = y_sum(1:dt:T, 1:dt:Lx);  colormap(hsv);
mesh(x, y, z, 'LineWidth', lineWidth);
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('$\gamma_{\hat{m},\hat{n}}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
grid on;  grid minor;  box on;  axis([0 T 0 Lx 0.3 0.38]);
%% u1 轨迹演化图
% figure('name', 'u1 轨迹演化图');
% [x,y] = meshgrid(1:dt:T, 1:dx:Lx);  z1 = u1_array(1:dt:T, 1:dt:Lx);  colormap(hsv);
% mesh(x, y, z1, 'LineWidth', lineWidth);
% set(gca,'FontSize', fontSizeAxis, 'Linewidth', lineWidth, 'FontWeight', 'bold');
% xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
% ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
% grid on;  grid minor;  box on;
% %% u2 轨迹演化图
% figure('name', 'u2 轨迹演化图');
% [x,y] = meshgrid(1:dt:T, 1:dx:Lx);  z2 = u2_array(1:dt:T, 1:dt:Lx);  colormap(hsv);
% mesh(x, y, z2, 'LineWidth', lineWidth);
% set(gca,'FontSize', fontSizeAxis, 'Linewidth', lineWidth, 'FontWeight', 'bold');
% xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
% ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
% grid on;  grid minor;  box on;
% %% u3 轨迹演化图
% figure('name', 'u3 轨迹演化图');
% [x,y] = meshgrid(1:dt:T, 1:dx:Lx);  z3 = u3_array(1:dt:T, 1:dt:Lx);  colormap(hsv);
% mesh(x, y, z3, 'LineWidth', lineWidth);
% set(gca,'FontSize', fontSizeAxis, 'Linewidth', lineWidth, 'FontWeight', 'bold');
% xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
% ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
% grid on;  grid minor;  box on;
% %% u4 轨迹演化图
% figure('name', 'u3 轨迹演化图');
% [x,y] = meshgrid(1:dt:T, 1:dx:Lx);  z4 = u4_array(1:dt:T, 1:dt:Lx);  colormap(hsv);
% mesh(x, y, z4, 'LineWidth', lineWidth);
% set(gca,'FontSize', fontSizeAxis, 'Linewidth', lineWidth, 'FontWeight', 'bold');
% xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
% ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeXY);
% grid on;  grid minor;  box on;
%% 传输信号减少率
disp('****************************************');
disp('*********H∞噪声抑制水平/触发参数*********');
disp('****************************************');
Gamma = gamma,  Trigger_thresholds = beta
disp('****************************************');
disp('************水平/垂直 触发次数***********');
disp('****************************************');
triggerTimes_h,  triggerTimes_v
disp('****************************************');
disp('*************水平/垂直 传输率************');
disp('****************************************');
Tans_rate_h = triggerTimes_h / (Nx*Nt)
Tans_rate_v = triggerTimes_v / (Nx*Nt)
disp('****************************************');
disp('************水平/垂直 平均传输率**********');
disp('****************************************');
Mean_Tans_rate_h = (sumTriggerTimes_h / numOfRepeat) / (Nx*Nt)
Mean_Tans_rate_v = (sumTriggerTimes_v / numOfRepeat) / (Nx*Nt)
disp('****************************************');
disp('***************爱趣无穷范数**************');
disp('****************************************');
Norm_y = sqrt(norm(y_array(1:Nt,1:Nx),2)) / sqrt(norm(w_array(1:Nt,1:Nx), 2))
disp(['运行时间: ', num2str(toc)]); 