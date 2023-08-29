function [x1,y1,x2,y2]=extraction()
global N_fig 
N_fig=3000;
r11=read_wobj('40.obj');%读取obj中数据
r22=read_wobj('60.obj');%读取obj中数据
%轮廓曲线1
coef1=r11.vertices;
coefs1=coef1';
coefs1(3,:)=[];
knots1=r11.objects(4).data;
leng1=length(knots1);
knots_l=ones(1,leng1);
for i=1:leng1
    knots_l(i)=(knots1(i)-knots1(1))/(knots1(leng1)-knots1(1));
end
line1 = nrbmak(coefs1,knots_l); 
ut1=linspace(0,1,N_fig);%在0~1间等分2000份
p1=nrbeval(line1,ut1); %带入求解公式得到对应点的坐标p
x1=p1(1,:);
y1=abs(p1(2,:));
%画图
% nrbplot(line1,2000);
% hold on
% plot(p1(1,:),p1(2,:),'--');
%修整坐标
x1=0.001*x1; %单位变换
y1=0.001*y1; %单位变换


%过渡曲线2
coef2=r22.vertices;
coefs2=coef2';
coefs2(3,:)=[];
knots2=r22.objects(4).data;
leng2=length(knots2);
knots_2=ones(1,leng2);
for i=1:leng2
    knots_2(i)=(knots2(i)-knots2(1))/(knots2(leng2)-knots2(1));
end
line2 = nrbmak(coefs2,knots_2); 
ut2=linspace(0,1,N_fig);%在0~1间等分4000份
p2=nrbeval(line2,ut2); %带入求解公式得到对应点的坐标p
x2=fliplr(p2(1,:));
y2=fliplr(abs(p2(2,:)));
%画图
% nrbplot(line2,4000);
% hold on
% plot(p2(1,:),p2(2,:),'--');
%修整坐标
x2=0.001*x2; %单位变换
y2=0.001*y2; %单位变换
