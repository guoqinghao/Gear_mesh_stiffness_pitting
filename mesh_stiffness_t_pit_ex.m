function [K_a,K_b,K_s]=mesh_stiffness_t_pit_ex(r_b,teta_b,alfa_1,alfa_2,leng,tc,cell_pit)
global E G L
global x11 y11
%初始化刚度，设置点的数量
K_a=zeros(1,leng);
K_b=zeros(1,leng);
K_s=zeros(1,leng);

if tc==1
    xdbeta=(alfa_2-alfa_1)/leng;%将啮合区间分成leng份，求出每份大小
%     A_x_test=cell(1,leng);
%     d_A_x_test=cell(1,leng);
%     I_x_test=cell(1,leng);
    for i=1:leng
        % 赋初值
        invK_a=0;
        invK_b=0;
        invK_s=0;
        beta=(i-1)*xdbeta+alfa_1; % 当前微分点角度
        xx=r_b*((beta+teta_b)*sin(beta)+cos(beta));
        [~,Index] = min(abs(x11-xx));%Index返回索引
        d=x11(Index); % 当前x坐标
        h=y11(Index); % 当前y坐标
        
        for j=1:Index-1  % 依次遍历每个点
            x1=x11(j);
            h1=y11(j); 
            x2=x11(j+1);
            h2=y11(j+1); 
            dx=abs(x2-x1);
            if isempty(cell_pit{2,j})  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 健康部分 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               

                A_x=(h1+h2)*L;
                I_x=1/12*(h1+h2)^3*L;                      
            else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 点蚀部分 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                cell_pit_now=cell_pit{2,j};
                % 求截面积
                d_A_x=sum(cell_pit_now(1,:).*cell_pit_now(2,:));
                yy=(h1+h2)/2;
                A_x=2*L*yy-d_A_x;
                
                % 求新的质心轴
                d_moment1=(L-cell_pit_now(1,:)).*cell_pit_now(2,:).*cell_pit_now(3,:);
                moment1=sum(d_moment1);
                N_moment=200; % 确定份数
                moment3=0;
                for n=1:N_moment
                    moment3=moment3+L*(yy-cell_pit_now(3,end)+...
                        cell_pit_now(2,end)/2)/N_moment*...
                       (cell_pit_now(3,end)+n*(yy-cell_pit_now(3,end)+...
                        cell_pit_now(2,end)/2)/N_moment);
                end                
                moment=moment1-moment3;
                ddx=abs(moment/A_x); % 新质心轴距中心轴距离

                % 求惯性矩
                d_moi1=(L-cell_pit_now(1,:)).*cell_pit_now(2,:).*...
                    (cell_pit_now(3,:)+ddx).^2;
                moi1=sum(d_moi1);
                N_moi=200; % 确定份数
                moi2=0;
                moi3=0;
                for n=1:N_moi
                    moi2=moi2+L*(cell_pit_now(3,end)+ddx)/N_moi*...
                        (n*(cell_pit_now(3,end)+ddx)/N_moi-(cell_pit_now(3,end)+ddx)/(2*N_moi))^2;
                end
                for n=1:N_moi
                    moi3=moi3+L*(yy-ddx)/N_moi*(n*(yy-ddx)/N_moi-(yy-ddx)/(2*N_moi))^2;
                end   
                I_x=moi1+moi2+moi3;
            

            end

            invK_a=invK_a+ (sin(beta)^2/(E*A_x))*dx;
            invK_b=invK_b+ (((d-x1)*cos(beta)-h*sin(beta))^2/(E*I_x))*dx;
            invK_s=invK_s+ (1.2*cos(beta)^2/(G*A_x))*dx;
        
        end
        
        K_a(i)=1/invK_a;
        K_b(i)=1/invK_b;
        K_s(i)=1/invK_s;
    end
end
end

