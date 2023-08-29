function [get_p1,dd,if_s]=extraction_pinion_f(stage,num)
%% 曲面
N_plane=400;
alfa= 1;    %密度系数
get_p1={};
dd=[];
if_s=ones(1,num);% 点蚀是否重叠
is_d=[];    % 重叠点蚀索引
if_s(is_d)=0;
% cc=30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%读取点蚀%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nu=stage*100+1:stage*100+num
    rr=read_wobj([num2str(nu) '.obj']);%读取obj中数据
    coef0=rr.vertices;
    %求cc
    coef0cc=coef0(:,1);
    diff_coef0cc=diff(coef0cc);
    for i=1:length(diff_coef0cc)
        if diff_coef0cc(i)<0
            cc=i;
            break
        end
    end

    x_cof0=coef0(:,1);
    y_cof0=coef0(:,3);
    z_cof0=coef0(:,2);
    coef0=[x_cof0';-y_cof0';z_cof0'];
    cc_dim=length(coef0(1,:))/cc;
    coefs0=ones(3,cc,cc_dim);
    for k=1:cc_dim
        coefs0(:,:,k)=coef0(:,(k+(cc-1)*(k-1)):(k+(cc-1)*(k-1)+(cc-1)));
    end
    coefs0=0.001.*coefs0;

    knotsu=rr.objects(4).data;
    lengu=length(knotsu);
    knots_u=ones(1,lengu);
    for i=1:lengu
        knots_u(i)=(knotsu(i)-knotsu(1))/(knotsu(lengu)-knotsu(1));
    end
    knotsv=rr.objects(5).data;
    lengv=length(knotsv);
    knots_v=ones(1,lengv);
    for i=1:lengv
        knots_v(i)=(knotsv(i)-knotsv(1))/(knotsv(lengv)-knotsv(1));
    end
    knots0{1}=knots_u;
    knots0{2}=knots_v;

    plane = nrbmak(coefs0,knots0); 
    ut1=linspace(0,1,N_plane);%在0~1间等分N_plane份
    vt1=linspace(0,1,N_plane);%在0~1间等分N_plane份
    p1=nrbeval(plane,{ut1 vt1});%带入求解公式得到对应点的坐标p
    get_p1{end+1}=[p1(1,:);p1(2,:);p1(3,:)];

    % 计算各点之间的距离
    % 横向距离（第三维）
    horizontal_delta=diff(p1,1,3); % 沿着第3维做1次差分
    horizontal_delta_sqr=sqrt(horizontal_delta(1,:).^2+horizontal_delta(2,:).^2+...
        horizontal_delta(3,:).^2);
    % 纵向距离（第二维）
    portrait_delta=diff(p1,1,2); % 沿着第2维做1次差分
    portrait_delta_sqr=sqrt(portrait_delta(1,:).^2+portrait_delta(1,:).^2+...
        portrait_delta(1,:).^2);
    delta_sqr=sqrt((mean(horizontal_delta_sqr)).^2+(mean(portrait_delta_sqr)).^2);
    dd(end+1)=alfa*delta_sqr;

    % %单位变换
    % x_p1=0.001*x_p1;
    % y_p1=0.001*y_p1;
    % z_p1=0.001*z_p1;
%     画图
%     nrbplot(plane,[20 20]);
%     hold on
end
