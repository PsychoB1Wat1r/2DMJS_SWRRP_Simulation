%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% 开环系统绘图代码 + 2024-4-3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%% 2D MJSs + 滑动窗口轮询协议（无零阶保持器） + 异步 %%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
heat_exchangers_jump % 导入系统参数
%% 绘图坐标 绘图线条 坐标轴字体 大小
fontSizeXY = 10;  lineWidth = 1;  fontSizeAxis = 15; markerSize = 9;  rand('state', 50);  randn('state', 50);
%% 绘图参数
setp = 40;  % 步长
T = setp;  dt = 1;  Nt = T/dt;    % 参数：水平方向
Lx = setp;  dx = 1;  Nx = Lx/dx;  % 参数：垂直方向
xh_array = zeros(Nt, Nx); % xh系统状态信号存储数组
xv_array = zeros(Nt, Nx); % xv系统状态信号存储数组
w_array = zeros(Nt, Nx);  % 噪声信号存储数组
%% 系统边界条件  2019 Neurocomputing Dehao Li Jinling Liang Fan Wang
xh_array(1:Nt, 1) = 4.2*cos(2:Nt+1)'*sin(1);
xv_array(1:Nt, 1) = 4.2*cos(2:Nt+1)'*sin(2);
xh_array(1, 1:Nx) = 4.2*tanh(1)*sin(2:Nx+1)';
xv_array(1, 1:Nx) = 4.2*tanh(2)*sin(2:Nx+1)';
%% 噪声信号
for i = 1 : Nt
    for j = 1 : Nx
         w_array(i, j) = -7/((i+1)*(j+1));
    end
end
%% Roesser-type 模型双向切换序列生成
SysSwitSeq = ones(Nt+1, Nx+1);   % 分段齐次马氏链 高级转移概率 模式序列存储矩阵
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
ConSwitSeq = ones(Nt+1, Nx+1); % 对偶隐马尔科夫模型 异步行为 1 模式序列存储矩阵
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
%% 绘图
for i = 1 : Nt
    for j = 1 : Nx
        if SysSwitSeq(i,j) == 1 % 系统模式切换
            As = A1;  Es = E1; 
        else
            As = A2;  Es = E2; 
        end
        % 系统的动态方程
        xh_array(i+1,j) = As(1,1)*xh_array(i,j) + As(1,2)*xv_array(i,j) + Es(1,1)*w_array(i,j);
        xv_array(i,j+1) = As(2,1)*xh_array(i,j) + As(2,2)*xv_array(i,j) + Es(2,1)*w_array(i,j);
    end
end
% load('xh_ol.mat'); load('xv_ol.mat'); % 导入数据
%% 水平方向系统开环演化图
figure('name', '水平方向系统开环演化图')
[x,y] = meshgrid(1:dt:T, 1:dx:Lx);
z = xh_array(1:Nt, 1:Nx);  surf(x, y, z, 'LineWidth', lineWidth, 'FaceAlpha', 1);  
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('$x^{h}_{m,n}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
axis([0 T 0 Lx -40 40]);  grid on;  grid minor;  box on;  colormap(copper);
% view(-24.7,31.8);
%% 垂直方向系统开环演化图
figure('name', '垂直方向系统开环演化图')
[x,y] = meshgrid(1:dt:T, 1:dx:Lx);
z = xv_array(1:Nt, 1:Nx);  surf(x, y, z, 'LineWidth', lineWidth, 'FaceAlpha', 1);  
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$m$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('$x^{v}_{m,n}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
axis([0 T 0 Lx -40 40]);  grid on;  grid minor;  box on;  colormap(copper);
% view(-24.7,31.8);
disp(['运行时间: ',num2str(toc)]); 