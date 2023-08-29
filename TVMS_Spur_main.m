tic
%% 清空工作区和命令窗口
clc
clear
close all

%% 输入齿轮参数
global T0 eps_alfa F delta_NSTE r1 r2
[r1,r2,r_a1,r_a2,r_f1,r_f2,r_b1,r_b2,teta_b1,teta_b2,...
    hf1,hf2,Le1,Le2,r_d1,r_d2,ep1,ep2,eg1,eg2]=parameter_setting();


%% 提取齿轮轮廓及点蚀坐标
global x11 y11 x22 y22
global N_fig %齿廓轮廓点

% 提取齿廓点云
[x11,y11,x22,y22]=extraction();

% 轻度点蚀，提取点蚀点云
% stage=1; % 第几阶段
% num=11; % 点蚀数量
% [get_p1,dd1,if_s1]=extraction_pinion_f(stage,num);
% get_p1为各个图形点云坐标数组,dd为点云最大相邻距离(分辨率)

% % 中度点蚀
% stage=2; % 第几阶段
% num=47; % 点蚀数量
% [get_p1,dd1,if_s1]=extraction_pinion_f(stage,num);

% % 重度点蚀
stage=3; % 第几阶段
num=81; % 点蚀数量
[get_p1,dd1,if_s1]=extraction_pinion_f(stage,num);


%% 点云切片及切条
% 创建轮廓切片数组
cell_pit1=cell(3,N_fig); % 相应位置切片点云和切条信息
% cell_pit2=cell(3,N_fig); % 相应位置切片点云和切条信息
% cell_pit3=cell(3,N_fig); % 相应位置切片点云和切条信息

% 轻度点蚀，切片处理
[cell_pit1]=slice_cutting(get_p1,if_s1,cell_pit1);
% % 中度点蚀，切片处理
% [cell_pit2]=slice_cutting(get_p2,if_s2,cell_pit2);
% % 重度点蚀，切片处理
% [cell_pit3]=slice_cutting(get_p3,if_s3,cell_pit3);

% 轻度点蚀，切条处理
[cell_pit1,get_strip_point1]=strip_cutting(cell_pit1);
% % 中度点蚀，切条处理
% [cell_pit2,get_strip_point2]=strip_cutting(cell_pit2);
% % 重度点蚀，切条处理
% [cell_pit3,get_strip_point3]=strip_cutting(cell_pit3);


%% 计算啮合极限坐标
[alfa_01,alfa_11,alfa_02,alfa_12,xx_01,xx_11,xx_02,xx_12,yy_01,yy_11,...
    yy_02,yy_12,alfag1_double1,alfag1_double2,alfap1_double1,alfap1_double2,...
    xg1_double1,xg1_double2,xp1_double1,xp1_double2]=...
    critical_point(r_b1,r_b2,teta_b1,teta_b2,r_a1,r_a2,Le1,Le2);


%% 计算啮合刚度分量
[teta_f1_ex,S_f1_ex,teta_f2_ex,S_f2_ex,length1,length2,length3,lens,k,...
          K_t1,K_t2,K_t3,K_f1,K_f2,K_f3,K_f_12,K_f_21,C_h1,C_h2,C_h3]=...
            mesh_stiffness(r_f1,r_f2,r_b1,teta_b1,alfa_01,alfa_11,r_b2,teta_b2,alfa_02,alfa_12,...
            hf1,alfag1_double1,hf2,alfap1_double1,alfag1_double2,alfap1_double2,xx_11,xx_01,cell_pit1);
        
% [teta_f1_ex,S_f1_ex,teta_f2_ex,S_f2_ex,length1,length2,length3,lens,k,...
%           K_t1,K_t2,K_t3,K_f1,K_f2,K_f3,K_f_12,K_f_21,C_h1,C_h2,C_h3]=...
%             mesh_stiffness(r_f1,r_f2,r_b1,teta_b1,alfa_01,alfa_11,r_b2,teta_b2,alfa_02,alfa_12,...
%             hf1,alfag1_double1,hf2,alfap1_double1,alfag1_double2,alfap1_double2,xx_11,xx_01,cell_pit);

       
%% 计算非线性力和传递误差
% m代表p的点蚀程度 1-健康,3-微弱点蚀,5-中度点蚀,7-严重点蚀
% n代表g的点蚀程度 2-健康,4-微弱点蚀,6-中度点蚀,8-严重点蚀
% p代表点蚀程度 1-单齿健康,2-单齿微弱点蚀
%% 计算负载的分担比
% 健康
mm1=1; % 齿对1
nn1=2; % 齿对2
mm2=1; % 齿对1
nn2=2; % 齿对2
pp1=1;
pp2=1;
[F1,F2,delta]=newton_method(K_t1,K_f1,K_f_12,C_h1,K_t3,K_f3,K_f_21,C_h3,length1,k,ep1,eg1,ep2,eg2,T0,r_b1,mm1,nn1,mm2,nn2,pp1,pp2);
F11(1,:)=F1;
F12(1,:)=F2;
delta0(1,:)=delta;

% 微弱点蚀1
mm1=3; % 齿对1
nn1=2; % 齿对2
mm2=1; % 齿对1
nn2=2; % 齿对2
pp1=2;
pp2=1;
[F1,F2,delta]=newton_method(K_t1,K_f1,K_f_12,C_h1,K_t3,K_f3,K_f_21,C_h3,length1,k,ep1,eg1,ep2,eg2,T0,r_b1,mm1,nn1,mm2,nn2,pp1,pp2);
F11(2,:)=F1;
F12(2,:)=F2;
delta0(2,:)=delta;
% 微弱点蚀2
mm1=1; % 齿对1
nn1=2; % 齿对2
mm2=1; % 齿对1
nn2=2; % 齿对2
pp1=1;
pp2=1;
[F1,F2,delta]=newton_method(K_t1,K_f1,K_f_12,C_h1,K_t3,K_f3,K_f_21,C_h3,length1,k,ep1,eg1,ep2,eg2,T0,r_b1,mm1,nn1,mm2,nn2,pp1,pp2);
F11(3,:)=F1;
F12(3,:)=F2;
delta0(3,:)=delta;


%% 计算总刚度
% 健康
%第一阶段
K1_last(1,:)=F./(delta0(1,:)-delta_NSTE); %此为同时考虑非线性接触刚度和齿间耦合对基体刚度的影响
%第二阶段
delta_h=C_h2(1,:)*F^k; 
K_H=F./delta_h;
K2_last(1,:)=1./(K_t2(1,:)+K_t2(2,:)+K_f2(1,:)+K_f2(2,:)+1./K_H);
%第三阶段
K3_last(1,:)=K1_last(1,:);
% 负载分担比
LSR1(1,:)=F11(1,:)./(F11(1,:)+F12(1,:));
LSR2(1,:)=F12(1,:)./(F11(1,:)+F12(1,:));
K_last(1,:)=[K1_last(1,:) K2_last(1,:) K3_last(1,:)];

% 微弱点蚀
%第一阶段
K1_last(2,:)=F./(delta0(2,:)-delta_NSTE); %此为同时考虑非线性接触刚度和齿间耦合对基体刚度的影响
%第二阶段
delta_h_pit1=C_h2(2,:)*F^k; 
K_H_pit1=F./delta_h_pit1;
K2_last(2,:)=1./(K_t2(3,:)+K_t2(2,:)+K_f2(1,:)+K_f2(2,:)+1./K_H_pit1);
%第三阶段
K3_last(2,:)=F./(delta0(3,:)-delta_NSTE); %此为同时考虑非线性接触刚度和齿间耦合对基体刚度的影响
%下面计算负载分担比
LSR1(2,:)=F11(2,:)./(F11(2,:)+F12(2,:));
LSR2(2,:)=F12(3,:)./(F11(3,:)+F12(3,:));
% 逆序fliplr
K_last(2,:)=[K1_last(2,:) K2_last(2,:) K3_last(2,:)];


%% 绘图
% dbeta=360*eps_alfa/(N1*(length2+length1-1)); %将齿轮作用角（齿轮从开始啮合到终止啮合转过的角度）分为(PTH-1)份
% PTH=lens;
% for i=1:PTH
%     betai(i)=dbeta*(i-1);
% end

xx=linspace(0,eps_alfa,lens);


% 总啮合刚度
xxlabel='meshing cycle';
yylabel='mesh stiffness(\itN/mm)';
plot_KF(xx,K_last,xxlabel,yylabel);

% kt啮合刚度
K_single(1,:)=1./([K_t1(1,:),K_t2(1,:),K_t3(1,:)]+[K_t1(2,:),K_t2(2,:),K_t3(2,:)]);
K_single(2,:)=1./([K_t1(3,:),K_t2(3,:),K_t3(3,:)]+[K_t1(2,:),K_t2(2,:),K_t3(2,:)]);
xxlabel='meshing cycle';
yylabel='mesh stiffness(\itN/mm)';
plot_KF(xx,K_single,xxlabel,yylabel);

% 负载分担比
LSRgp(1,:)=[LSR1(1,:) ones(1,length2) LSR2(1,:)];
LSRgp(2,:)=[LSR1(2,:) ones(1,length2) LSR2(2,:)];
xxlabel='meshing cycle';
yylabel='Load sharing ratio';
plot_KF(xx,LSRgp,xxlabel,yylabel);

% % F11负载变化
% xxx=linspace(0,eps_alfa,length1);
% xxlabel='meshing cycle';
% yylabel='F11';
% plot_FF(xxx,F11,xxlabel,yylabel);
% % F12负载变化
% xxlabel='meshing cycle';
% yylabel='F12';
% plot_FF(xxx,F12,xxlabel,yylabel);


toc
% H=K_last(1,:);SL=K_last(2,:);
% save('Kt_pit.mat','H','SL')
% MD=K_last(2,:);
% save('Kt_pit.mat', 'MD', '-append')