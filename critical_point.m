function [alfa_01,alfa_11,alfa_02,alfa_12,xx_01,xx_11,xx_02,xx_12,yy_01,yy_11,...
    yy_02,yy_12,alfag1_double1,alfag1_double2,alfap1_double1,alfap1_double2,...
    xg1_double1,xg1_double2,xp1_double1,xp1_double2]=...
    critical_point(r_b1,r_b2,teta_b1,teta_b2,r_a1,r_a2,Le1,Le2)

global N1 N2 r1 r2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算极限角度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms alfa_01 alfa_11 alfa_02 alfa_12
alfa_01=bisection_method...
    (r_b1*(alfa_01+teta_b1)-Le1,-teta_b1,pi/2,1e-10);                %齿轮1起始极限点角度
alfa_11=bisection_method...
	(r_b1*(alfa_11+teta_b1)-sqrt(r_a1^2-r_b1^2),-teta_b1,pi/2,1e-10);%齿轮1分离极限点角度

alfa_02=bisection_method...
    (r_b2*(alfa_02+teta_b2)-sqrt(r_a2^2-r_b2^2),-teta_b2,pi/2,1e-10);%齿轮2起始极限点角度
alfa_12=bisection_method...
    (r_b2*(alfa_12+teta_b2)-Le2,-teta_b2,pi/2,1e-10);                %齿轮2分离极限点角度


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算对应的极限坐标%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%齿轮1
xx_01=r_b1*cos(alfa_01)+Le1*sin(alfa_01);
yy_01=Le1*cos(alfa_01)-r_b1*sin(alfa_01);
xx_11=r_b1*cos(alfa_11)+sqrt(r_a1^2-r_b1^2)*sin(alfa_11);
yy_11=sqrt(r_a1^2-r_b1^2)*cos(alfa_11)-r_b1*sin(alfa_11);

%齿轮2
xx_02=r_b2*cos(alfa_02)+sqrt(r_a2^2-r_b2^2)*sin(alfa_02);
yy_02=sqrt(r_a2^2-r_b2^2)*cos(alfa_02)-r_b2*sin(alfa_02);
xx_12=r_b2*cos(alfa_12)+Le2*sin(alfa_12);
yy_12=Le2*cos(alfa_12)-r_b2*sin(alfa_12);


%%%%%%%%%%%%%%%%%%%%%%%%计算单双齿啮合区间对应角度%%%%%%%%%%%%%%%%%%%%%%%%%%
syms alfag1_double1
% 齿轮1-齿对1
% alfa_01---alfag1_double1为双齿啮合阶段1
% alfag1_double1---alfag1_double2为单齿啮合阶段
% alfag1_double2---alfa_11为双齿啮合阶段2
alfag1_double1=bisection_method...
    (r_b1*(alfag1_double1+2*pi/N1+teta_b1)-sqrt(r_a1^2-r_b1^2),-teta_b1,pi/2,1e-10); % 齿对1啮合阶段临界点1
deltag=alfag1_double1-alfa_01;
alfag1_double2=alfa_11-deltag;                   % 齿对1啮合阶段临界点2


% 齿轮2-齿对1
% alfa_02---alfap1_double1为双齿啮合阶段1
% alfap1_double1---alfap1_double2为单齿啮合阶段
% alfap1_double2---alfa_12为双齿啮合阶段2
deltap=deltag*N1/N2;
alfap1_double1=alfa_02-deltap;                   % 齿对2啮合阶段临界点1
alfap1_double2=alfa_12+deltap;                   % 齿对2啮合阶段临界点2



%%%%%%%%%%%%%%%%%%计算单双齿啮合区间对应坐标(未使用)%%%%%%%%%%%%%%%%%%%%%%%%
%齿轮1
xg1_double1=r_b1*((alfag1_double1+teta_b1)*sin(alfag1_double1)+cos(alfag1_double1));
xg1_double2=r_b1*((alfag1_double2+teta_b1)*sin(alfag1_double2)+cos(alfag1_double2));


%齿轮2
xp1_double1=r_b2*((alfap1_double1+teta_b2)*sin(alfap1_double1)+cos(alfap1_double1));
xp1_double2=r_b2*((alfap1_double2+teta_b2)*sin(alfap1_double2)+cos(alfap1_double2));
