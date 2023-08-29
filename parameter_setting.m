function [r1,r2,r_a1,r_a2,r_f1,r_f2,r_b1,r_b2,teta_b1,teta_b2,...
    hf1,hf2,Le1,Le2,r_d1,r_d2,ep1,ep2,eg1,eg2]=parameter_setting()
%% 输入齿轮参数
global N1 N2 m alfa h_a_star c_star r_int1 r_int2 E Nu L T0 G...
    eps_alfa F delta_NSTE inv_alfa
m=2e-3;                    % 模数，单位为m
N1=40;                     % 齿轮1齿数
N2=60;                     % 齿轮2齿数
alfa=20*2*pi/360;          % 压力角
h_a_star=1;                % 齿顶系数
c_star=0.25;               % 顶隙系数
E=2.05e11;                 % 弹性模量
Nu=0.3;                    % 泊松比
L=20e-3;                   % 齿宽，单位为m
r_int1=20e-3;              % 轮毂孔半径1
r_int2=30e-3;              % 轮毂孔半径2
T0=40;                     % 输入转矩Nm
ep1=0.00e-3;               % 齿对1驱动侧齿廓误差
eg1=0.00e-3;               % 齿对1从动侧齿廓误差
ep2=0.00e-3;               % 齿对2驱动侧齿廓误差
eg2=0.00e-3;               % 齿对2从动侧齿廓误差

%% 计算初始量
G=E/(2*(1+Nu));                         % 切变模量
r1=m*N1/2;                              % 分度圆半径1
r_a1=r1+h_a_star*m;                     % 齿顶圆半径1
r_f1=r1-(c_star+h_a_star)*m;            % 齿根圆半径1
r_b1=r1*cos(alfa);                      % 基圆半径
r2=m*N2/2;                              % 分度圆半径2
r_a2=r2+h_a_star*m;                     % 齿顶圆半径2
r_f2=r2-(c_star+h_a_star)*m;            % 齿根圆半径2
r_b2=r2*cos(alfa);                      % 基圆半径2
eps_alfa=(sqrt(r_a2^2-r_b2^2)+sqrt(r_a1^2-r_b1^2)-(r1+r2)*sin(alfa))/(pi*m*cos(alfa));
                                        % 重合度计算，B1B2/Pb
inv_alfa=tan(alfa)-alfa;                % 渐开线函数
teta_b1=pi/(2*N1)+inv_alfa;             % 基圆上半齿角1
teta_b2=pi/(2*N2)+inv_alfa;             % 基圆上半齿角2

% 用于计算基体刚度
hf1=r_f1/r_int1;                        % 齿根圆半径1/轮毂孔半径1
hf2=r_f2/r_int2;                        % 齿根圆半径2/轮毂孔半径2
 
% 起始啮合圆
Le1=sqrt((r1+r2)^2-(r_b1+r_b2)^2)-sqrt(r_a2^2-r_b2^2);
r_d1=sqrt(Le1^2+r_b1^2);
Le2=sqrt((r1+r2)^2-(r_b1+r_b2)^2)-sqrt(r_a1^2-r_b1^2);
r_d2=sqrt(Le2^2+r_b2^2);

F=T0/r_b1;                             % 啮合线上总啮合力
delta_NSTE=min(ep1+eg1,ep2+eg2);       % 无载荷静态传递误差