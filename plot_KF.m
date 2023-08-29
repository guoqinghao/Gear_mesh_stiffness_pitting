function []=plot_KF(xx,K_last,xxlabel,yylabel)
figure
p1=plot(xx,K_last(1,:),'color',[0 0.2844 0.5111],'LineWidth',1.5);
hold on
plot(xx,K_last(2,:),'--','LineWidth',1.5);
hold on
g = get(p1,'Parent');%对应p1所在的坐标轴
set(g,'Linewidth',1.5,'FontSize',12,'FontName','Times New Roman','FontWeight','bold');
%这里定义坐标轴的线宽，字体，字号，字体是否加粗
ylabel(yylabel,'FontSize',12,'FontName','Times New Roman','FontWeight','bold');
xlabel(xxlabel,'FontSize',12,'FontName','Times New Roman','FontWeight','bold');
%\it表示斜体，\rm表示正体