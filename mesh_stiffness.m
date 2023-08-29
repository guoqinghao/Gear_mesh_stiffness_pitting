function [teta_f1_ex,S_f1_ex,teta_f2_ex,S_f2_ex,length1,length2,length3,lens,k,...
          K_t1,K_t2,K_t3,K_f1,K_f2,K_f3,K_f_12,K_f_21,C_h1,C_h2,C_h3]=...
            mesh_stiffness(r_f1,r_f2,r_b1,teta_b1,alfa_01,alfa_11,r_b2,teta_b2,alfa_02,alfa_12,...
            hf1,alfag1_double1,hf2,alfap1_double1,alfag1_double2,alfap1_double2,xx_11,xx_01,cell_pit)
global x11 y11 y22 E L eps_alfa N1 N2  
%% 计算齿根圆半齿角及齿厚
%齿轮1
    teta_f1_ex=asin(y11(1)/r_f1);   % 齿根半齿角1
    S_f1_ex=2*y11(1);               % 齿根圆对应齿厚1

%齿轮2
    teta_f2_ex=asin(y22(1)/r_f2);   % 齿根半齿角2
    S_f2_ex=2*y22(1);               % 齿根圆对应齿厚2

    
%% 计算啮合刚度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%分点%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
len=200; %设置单周期总长度
length1=floor(len*(eps_alfa-1)/eps_alfa); %第一阶段需要的最终长度
length2=floor(len*(1-(eps_alfa-1))/eps_alfa); %第二阶段需要的最终长度
length3=length1; %第三阶段需要的最终长度
lens=length1+length2+length3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%计算齿轮1的K_a,K_b,K_s%%%%%%%%%%%%%%%%%%%%%%%%%%
% K_t1,K_t2,K_t3分别代表阶段1，阶段2，阶段3
% K_t1=[g1;p1;g_pit11;p_pit11;g_pit21;p_pit21;g_pit31;p_pit31];
% K_t2=[g2;p2;g_pit12;p_pit12;g_pit22;p_pit22;g_pit32;p_pit32];
% K_t3=[g3;p3;g_pit13;p_pit13;g_pit23;p_pit23;g_pit33;p_pit33];
%健康 
[K_a,K_b,K_s]=mesh_stiffness_t_ex...
    (r_b1,teta_b1,alfa_01,alfa_11,lens,1); %单齿啮合区间
K_t=1./K_a+1./K_b+1./K_s;
K_t1(1,:) = K_t(1,1:length1);
K_t2(1,:) = K_t(1,(length1+1):(length2+length1));
K_t3(1,:) = K_t(1,(lens-length3+1):lens);

%点蚀1
[K_a,K_b,K_s]=mesh_stiffness_t_pit_ex...
    (r_b1,teta_b1,alfa_01,alfa_11,lens,1,cell_pit); %单齿啮合区间
K_t=1./K_a+1./K_b+1./K_s;
K_t1(3,:) = K_t(1,1:length1);
K_t2(3,:) = K_t(1,(length1+1):(length2+length1));
K_t3(3,:) = K_t(1,(lens-length3+1):lens);

%%%%%%%%%%%%%%%%%%%%%%%%%%%计算齿轮2的K_a,K_b,K_s%%%%%%%%%%%%%%%%%%%%%%%%%%
[K_a,K_b,K_s]=mesh_stiffness_t_ex...
    (r_b2,teta_b2,alfa_02,alfa_12,lens,2); %单齿啮合区间
K_t=1./K_a+1./K_b+1./K_s;
K_t1(2,:) = K_t(1,1:length1);
K_t2(2,:) = K_t(1,(length1+1):(length2+length1));
K_t3(2,:) = K_t(1,(lens-length3+1):lens);

%%%%%%%%%%%%%%%%%%计算齿轮1的K_f,K_f12,K_f21,K_f10,K_f01%%%%%%%%%%%%%%%%%%%
[K_f,K_f12,K_f21]=mesh_stiffness_f(r_b1,r_f1,teta_f1_ex,S_f1_ex,hf1,...
    teta_b1,alfa_01,alfag1_double1,N1,1,12,length1); %齿轮1齿对1啮合阶段1
K_f1(1,:)=1./K_f;
K_f_12(1,:)=1./K_f12;
K_f_21(1,:)=1./K_f21;

[K_f,~,~]=mesh_stiffness_f(r_b1,r_f1,teta_f1_ex,S_f1_ex,hf1,...
    teta_b1,alfag1_double1,alfag1_double2,N1,1,0,length2); %齿轮1齿对1啮合阶段2
K_f2(1,:)=1./K_f;

[K_f,~,~]=mesh_stiffness_f(r_b1,r_f1,teta_f1_ex,S_f1_ex,hf1,...
    teta_b1,alfag1_double2,alfa_11,N1,1,0,length3); %齿轮1齿对1啮合阶段2
K_f3(1,:)=1./K_f;

%%%%%%%%%%%%%%%%%%计算齿轮2的K_f,K_f12,K_f21,K_f10,K_f01%%%%%%%%%%%%%%%%%%%
[K_f,K_f12,K_f21]=mesh_stiffness_f(r_b2,r_f2,teta_f2_ex,S_f2_ex,hf2,...
    teta_b2,alfa_02,alfap1_double1,N2,2,12,length1); %齿轮2齿对1啮合阶段1
K_f1(2,:)=1./K_f;
K_f_12(2,:)=1./K_f12;
K_f_21(2,:)=1./K_f21;

[K_f,~,~]=mesh_stiffness_f(r_b2,r_f2,teta_f2_ex,S_f2_ex,hf2,...
    teta_b2,alfap1_double1,alfap1_double2,N2,2,0,length2); %齿轮2齿对1啮合阶段2
K_f2(2,:)=1./K_f;

[K_f,~,~]=mesh_stiffness_f(r_b2,r_f2,teta_f2_ex,S_f2_ex,hf2,...
    teta_b2,alfap1_double2,alfa_12,N2,2,0,length3); %齿轮2齿对1啮合阶段2
K_f3(2,:)=1./K_f;

% K_f_singleg1_12_ex=zeros(1,length1);
% K_f_singleg1_21_ex=zeros(1,length1);
% K_f_singlep1_12_ex=zeros(1,length1);
% K_f_singlep1_21_ex=zeros(1,length1);                                      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算齿轮1的C_h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 健康
% 健康
k=0.9;
C_h=1.275/(0.9.*E^0.9.*L^0.8);
% 区分单双齿
C_h1(1,:)=C_h.*ones(1,length1);
C_h2(1,:)=C_h.*ones(1,length2);
C_h3(1,:)=C_h.*ones(1,length3);

% 点蚀
% 获取整个齿廓的接触长度矩阵
cell_pit_L=zeros(1,lens);            % 创建接触长度矩阵
xdbeta=(alfa_11-alfa_01)/lens;                         %将啮合区间分成leng份，求出每份大小
for i=1:lens
    beta1=(i-1)*xdbeta+alfa_01;                        % 齿轮1当前微分点角度
    xx1=r_b1*((beta1+teta_b1)*sin(beta1)+cos(beta1));  % 齿轮1当前横坐标
    [~,Index1] = min(abs(x11-xx1));                    %Index1返回齿轮1当前横坐标对应索引
    L1=cell_pit{2,Index1};                            % 表面接触长度坐标1
    % 下面进行计算最终有效长度
    if ~isempty(L1) 
        cell_pit_L(i)= L1(1,1);
    end
end
L_now=L-cell_pit_L;






% cell_pit_L=zeros(1,length(cell_pit(1,:)));       % 创建接触长度矩阵
% Index_L=find( ~cellfun('isempty',cell_pit(2,:))); 
% for ii=1:length(Index_L)
%     i=Index_L(ii);
% %     if isempty(cell_pit{2,i})
% %         cell_pit_L(i)=0;
% %     else  
%         cell_pit_L(i)=cell_pit{2,i}(1,1);
% %     end
% end
% % 获取测量点的接触长度矩阵
% xdbeta=(xx_11-xx_01)/lens;%将啮合区间分成len份，求出每份大小
% Index_L2=zeros(1,lens);
% for i=1:lens
%     xbeta=(i-1)*xdbeta+xx_01; %当前微分点坐标
%     [~,Index_L1] = min(abs(x11-xbeta));%Index返回索引
%     Index_L2(i)=Index_L1;
% end
% % 获取测量点的C_h矩阵
% L_select=cell_pit_L(Index_L2);
% L_now=L-L_select;
C_h_pit=1.275./(0.9.*E.^0.9.*L_now.^0.8);
% plot(C_h_pit);
% 区分单双齿
C_h1(2,:)=C_h_pit(1:length1);
C_h2(2,:)=C_h_pit((length1+1):(length1+length2));
C_h3(2,:)=C_h_pit((length1+length2+1):lens);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算齿轮2的C_h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%健康
% C_h1(2,:)=C_h.*ones(1,length1);
% C_h2(2,:)=C_h.*ones(1,length2);
% C_h3(2,:)=C_h.*ones(1,length3);

