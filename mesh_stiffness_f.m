function [K_f,K_f12,K_f21]=mesh_stiffness_f(r_b,r_f,teta_f,S_f,hf,...
    teta_b,alfa_1,alfa_2,N,gear,tooth_pair,leng)
global E L
global x11 y11 x22 y22
%% 求齿轮基体柔性变形刚度用到的参数
Coef = [ -3.785e-5   1.931e-3  1.293e-3  -2.254e-2  9.922e-5   -1.306e-2  2.6e-4   0.2451  3.5260;
         -8.177e-5   4.312e-3  2.547e-3  3.4e-2     -2.184e-3  -4.419e-2  1.095e-3   0.6473  0.5298;
         -4.415e-5   1.669e-2  2.053e-3  -8.121e-3  -2.115e-4  6.338e-3   -1.318e-4  0.9658  0.4352;
         50.15       2.190     -4.765    -4.636     5.705      -7.071     11.26      3.4580  -0.1916;
         ];
Ai=Coef(:,1);    Bi=Coef(:,2);    Ci=Coef(:,3);
Di=Coef(:,4);    Ei=Coef(:,5);    Fi=Coef(:,6);
Gi=Coef(:,7);    Hi=Coef(:,8);    Ii=Coef(:,9);

for i=1:3
    Xi(i)=Ai(i)/teta_f^3 + Bi(i)*hf^3 + Ci(i)/teta_f^2 + Di(i)*hf^2 +...
          Ei(i)*hf^2/teta_f + Fi(i)*hf/teta_f + Gi(i)*hf/teta_f^2 +...
          Hi(i)*hf + Ii(i);
end 
L_star=Xi(1);    M_star=Xi(2);    P_star=Xi(3);    
Q_star=Ai(4)*teta_f^3 + Bi(4)/hf^3 + Ci(4)*teta_f^2 + Di(4)/hf^2 +...
       Ei(4)*teta_f/hf^2 + Fi(4)*teta_f/hf + Gi(4)*teta_f^2/hf +...
       Hi(4)/hf + Ii(4);

Coef1 = [-2.42e-5   1.607e-3  8.513e-4  -1.614e-2  -4.55e-5   -1.001e-2  2.076e-4   0.1931   -0.8056;
         40.21      -0.3763   -8.284    0.5267     -5.404     3.073      -12.77     -0.2074  0.5662;
         -3.620e-5  1.967e-3  1.074e-3  1.940e-2   -1.100e-3  -2.254e-2  5.622e-4   0.3177   -0.2317;
         -39.58     0.3751    8.117     -0.5241    5.376      -3.019     12.47      0.2050   -0.5785;
         -3.610e-5  1.960e-3  1.070e-3  1.952e-2   -1.101e-3  -2.252e-2  5.621e-4   0.3170   -0.2305;
         -3.570e-5  1.099e-2  1.689e-3  -0.1434    -3.802e-4  8.386e-3   -1.440e-4  0.6640   -0.7143;
         3.597e02   9.839e-2  -61.51e1  -0.1819    -2.452     0.8327     22.93      0.1655   0.1135;
         -3.564e02  -9.638e-2 60.86e1   0.1778     2.503      -0.8984    -22.53     -0.1627  -0.1146;
         -4.210e-5  1.638e-2  1.962e-3  -3.132e-3  -2.771e-4  7.486e-3   -1.463e-4  0.9336   -0.5324;
        ];
A1=Coef1(:,1);    B1=Coef1(:,2);    C1=Coef1(:,3);
D1=Coef1(:,4);    E1=Coef1(:,5);    F1=Coef1(:,6);
G1=Coef1(:,7);    H1=Coef1(:,8);    I1=Coef1(:,9);

for i=1:9
    X11(i)=A1(i)/teta_f^3 + B1(i)*hf^3 + C1(i)/teta_f^2 + D1(i)*hf^2 +...
           E1(i)*hf^2/teta_f + F1(i)*hf/teta_f + G1(i)*hf/teta_f^2 +...
           H1(i)*hf + I1(i);
end 
L1=X11(1);   P1=X11(3);   R1=X11(5);    S1=X11(6);  V1=X11(9);

for i=1:9
    X12(i)=A1(i)*teta_f^3 + B1(i)/hf^3 + C1(i)*teta_f^2 + D1(i)/hf^2 +...
           E1(i)*teta_f/hf^2 + F1(i)*teta_f/hf + G1(i)*teta_f^2/hf +...
           H1(i)/hf + I1(i);
end 
M1=X12(2);   Q1=X12(4);   T1=X12(7);    U1=X12(8);  

Coef2 = [-2.42e-5  1.607e-3   8.513e-4  -1.614e-2  -4.55e-5   -1.001e-2  2.076e-4  0.1931   -0.8056;
        -40.21     0.3763     8.284     -0.5267    5.404      -3.073     12.77     0.2074   -0.5662;
        -3.620e-5  1.967e-3   1.074e-3  1.940e-2   -1.100e-3  -2.254e-2  5.622e-4  0.3177   -0.2317;
        39.58      -0.3751    -8.117    0.5241     -5.376     3.019      -12.47    -0.2050  0.5785;
        -3.610e-5  1.960e-3   1.070e-3  1.952e-2   -1.101e-3  -2.252e-2  5.621e-4  0.3170   -0.2305;
        -3.570e-5  1.099e-2   1.689e-3  -0.1434    -3.802e-4  8.386e-3   -1.440e-4 0.6640   -0.7143;
        -3.597e02  -9.839e-2  61.51e1     0.1819   2.452      -0.8327    -22.93    -0.1655  -0.1135;
        3.564e02   9.638e-2   -60.86e1  -0.1778    -2.503     0.8984     22.53     0.1627   0.1146;
        -4.210e-5  1.638e-2  1.962e-3  -3.132e-3   -2.771e-4  7.486e-3   -1.463e-4 0.9336   -0.5324;
        ];
A2=Coef2(:,1);    B2=Coef2(:,2);    C2=Coef2(:,3);
D2=Coef2(:,4);    E2=Coef2(:,5);    F2=Coef2(:,6);
G2=Coef2(:,7);    H2=Coef2(:,8);    I2=Coef2(:,9);

for i=1:9
    X21(i)=A2(i)/teta_f^3 + B2(i)*hf^3 + C2(i)/teta_f^2 + D2(i)*hf^2 +...
           E2(i)*hf^2/teta_f + F2(i)*hf/teta_f + G2(i)*hf/teta_f^2 +...
           H2(i)*hf + I2(i);
end 
L2=X21(1);   P2=X21(3);   R2=X21(5);    S2=X21(6);  V2=X21(9);

for i=1:9
    X22(i)=A2(i)*teta_f^3 + B2(i)/hf^3 + C2(i)*teta_f^2 + D2(i)/hf^2 +...
        E2(i)*teta_f/hf^2 + F2(i)*teta_f/hf + G2(i)*teta_f^2/hf +...
        H2(i)/hf + I2(i);
end 
M2=X22(2);   Q2=X22(4);   T2=X22(7);    U2=X22(8); 


%% 计算基体刚度
K_f=zeros(1,leng);
K_f21=zeros(1,leng);
K_f12=zeros(1,leng);

xdbeta=(alfa_2-alfa_1)/leng;%将啮合区间分成leng份，求出每份大小
if gear==1
    for i=1:leng
        beta=(i-1)*xdbeta+alfa_1; % 当前微分点角度
        xx=r_b*((beta+teta_b)*sin(beta)+cos(beta));
        [~,Index] = min(abs(x11-xx));%Index返回索引
        d=x11(Index); % 当前x坐标
        h=y11(Index); % 当前y坐标
        u_f=d-r_f-h*tan(beta);
        %力作用线与中线交点到齿根距离（为了求基体刚度）
        invK_f=(cos(beta)^2/(E*L))*(L_star*(u_f/S_f)^2 +...
                M_star*(u_f/S_f) + P_star*(1+Q_star*(tan(beta))^2));
        K_f(i)=1/invK_f;    

        %计算基体变形耦合影响
        if tooth_pair==12 %判断齿对
            %齿牙1对齿牙2的基体变形耦合影响  
            if gear==1 %判断齿轮1和齿轮2，角度换算方法不同
                beta2=beta+2*pi/N; %当前微分点齿牙2角度
            else
                beta2=beta-2*pi/N;
            end
            d2=r_b*((beta2+teta_b)*sin(beta2)+cos(beta2)); %当前x坐标
            [~,Index2] = min(abs(x11-d2));%Index2返回索引
            h2=y11(Index2);
            u_f2=d2-r_f-h2*tan(beta2);
            %力作用线与中线交点到齿根距离（为了求基体刚度）
            invK_f21=cos(beta)*cos(beta2)*(L1*(u_f2*u_f/(S_f*S_f))^2+...
                (tan(beta2)*M1+P1)*u_f/S_f+(tan(beta)*Q1+R1)*u_f2/S_f+...
                (tan(beta)*S1+T1)*tan(beta2)+U1*tan(beta)+V1)/(E*L);
            invK_f12=cos(beta)*cos(beta2)*(L2*(u_f*u_f2/(S_f*S_f))^2+...
                (tan(beta)*M2+P2)*u_f2/S_f+(tan(beta2)*Q2+R2)*u_f/S_f+...
                (tan(beta2)*S2+T2)*tan(beta)+U2*tan(beta2)+V2)/(E*L);
            K_f21(i)=1/invK_f21;
            K_f12(i)=1/invK_f12;
        end

    end
else
    for i=1:leng
        beta=(i-1)*xdbeta+alfa_1; % 当前微分点角度
        xx=r_b*((beta+teta_b)*sin(beta)+cos(beta));
        [~,Index] = min(abs(x22-xx));%Index返回索引
        d=x22(Index); % 当前x坐标
        h=y22(Index); % 当前y坐标
        u_f=d-r_f-h*tan(beta);
        %力作用线与中线交点到齿根距离（为了求基体刚度）
        invK_f=(cos(beta)^2/(E*L))*(L_star*(u_f/S_f)^2 +...
                M_star*(u_f/S_f) + P_star*(1+Q_star*(tan(beta))^2));
        K_f(i)=1/invK_f;    

        %计算基体变形耦合影响
        if tooth_pair==12 %判断齿对
            %齿牙1对齿牙2的基体变形耦合影响  
            if gear==1 %判断齿轮1和齿轮2，角度换算方法不同
                beta2=beta+2*pi/N; %当前微分点齿牙2角度
            else
                beta2=beta-2*pi/N;
            end
            d2=r_b*((beta2+teta_b)*sin(beta2)+cos(beta2)); %当前x坐标
            [~,Index2] = min(abs(x22-d2));%Index2返回索引
            h2=y22(Index2);
            u_f2=d2-r_f-h2*tan(beta2);
            %力作用线与中线交点到齿根距离（为了求基体刚度）
            invK_f21=cos(beta)*cos(beta2)*(L1*(u_f2*u_f/(S_f*S_f))^2+...
                (tan(beta2)*M1+P1)*u_f/S_f+(tan(beta)*Q1+R1)*u_f2/S_f+...
                (tan(beta)*S1+T1)*tan(beta2)+U1*tan(beta)+V1)/(E*L);
            invK_f12=cos(beta)*cos(beta2)*(L2*(u_f*u_f2/(S_f*S_f))^2+...
                (tan(beta)*M2+P2)*u_f2/S_f+(tan(beta2)*Q2+R2)*u_f/S_f+...
                (tan(beta2)*S2+T2)*tan(beta)+U2*tan(beta2)+V2)/(E*L);
            K_f21(i)=1/invK_f21;
            K_f12(i)=1/invK_f12;
        end

    end    


end