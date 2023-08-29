function [K_a,K_b,K_s]=mesh_stiffness_t_ex(r_b,teta_b,alfa_1,alfa_2,leng,tc)
global E G L
global x11 y11 x22 y22
%初始化刚度，设置点的数量
K_a=zeros(1,leng);
K_b=zeros(1,leng);
K_s=zeros(1,leng);
%     A_x_testh=cell(1,leng);
%     I_x_testh=cell(1,leng);
if tc==1
    xdbeta=(alfa_2-alfa_1)/leng;%将啮合区间分成leng份，求出每份大小
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
        for k=1:(Index-1)
            x1=x11(k);
            h1=y11(k); 
            x2=x11(k+1);
            h2=y11(k+1); 
            dx=abs(x2-x1);
            A_x=(h1+h2)*L;
            I_x=1/12*(h1+h2)^3*L;
%             
%             A_x_testh{i}(end+1)=A_x;    
%             I_x_testh{i}(end+1)=I_x; 
            invK_a=invK_a+ (sin(beta)^2/(E*A_x))*dx;
            invK_b=invK_b+ (((d-x1)*cos(beta)-h*sin(beta))^2/(E*I_x))*dx;
            invK_s=invK_s+ (1.2*cos(beta)^2/(G*A_x))*dx;
        end
        K_a(i)=1/invK_a;
        K_b(i)=1/invK_b;
        K_s(i)=1/invK_s;   
    end
end

if tc==2
    xdbeta=(alfa_2-alfa_1)/leng;%将啮合区间分成leng份，求出每份大小
    for i=1:leng
        
        % 赋初值
        invK_a=0;
        invK_b=0;
        invK_s=0;
        beta=(i-1)*xdbeta+alfa_1; % 当前微分点角度
        xx=r_b*((beta+teta_b)*sin(beta)+cos(beta));
        [~,Index] = min(abs(x22-xx));%Index返回索引
        d=x22(Index); % 当前x坐标
        h=y22(Index); % 当前y坐标
        for k=1:(Index-1)
            x1=x22(k);
            h1=y22(k); 
            x2=x22(k+1);
            h2=y22(k+1); 
            dx=abs(x2-x1);
            A_x=(h1+h2)*L;
            I_x=1/12*(h1+h2)^3*L;
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