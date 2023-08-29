function [F1,F2,delta]=newton_method...
    (K_t1,K_f1,K_f_12,C_h1,K_t3,K_f3,K_f_21,C_h3,length1,k,ep1,eg1,ep2,eg2,T0,r_b1,mm1,nn1,mm2,nn2,pp1,pp2)
for i=1:length1  
    syms x y z                    % 定义符号变量x y z
    f1=x.*(K_t1(mm1,i)+K_t1(nn1,i)+K_f1(1,i)+K_f1(2,i))+...
        C_h1(pp1,i).*x.^k +y*(K_f_12(1,i)+K_f_12(2,i))+ep1+eg1-z; %方程1
    f2=y.*(K_t3(mm2,i)+K_t3(nn2,i)+K_f3(1,i)+K_f3(2,i))+...
        C_h3(pp2,i).*y.^k +x*(K_f_21(1,i)+K_f_21(2,i))+ep2+eg2-z; %方程2
    f3=x+y-T0./r_b1;              % 方程3
    f=[f1;f2;f3];                 % 三个方程组成方程组
    %%计算雅可比行列式
    J=[diff(f1,x) diff(f1,y) diff(f1,z);diff(f2,x) diff(f2,y) diff(f2,z);diff(f3,x) diff(f3,y) diff(f3,z);];
    %%自定义牛顿迭代法，求方程组的解
    n=1;                          % 记录迭代次数
    x0=[240;290;2e-5];            % 迭代初始值
    e=1;
    while e>1e-5                  % 如果精度不满足要求，则一直进行迭代，直到满足要求为止
        x=x0(1);                  % 给x赋初始值
        y=x0(2);                  % 给y赋初始值
        z=x0(3);                  % 给z赋初始值
        x1=x0-(eval(J))\eval(f);  % 迭代计算，用x0迭代得到x1
        e=max(abs((x1-x0)./x1));  % 计算求解精度
        x0=x1;                    % 将当前得到的x1赋值给x0，作为下一次迭代的初始值
        n=n+1;                    % 迭代次数加1
    end
    %%显示方程组的解，并将求解结果代回方程组，验算求解结果
    x=x1(1);                      % x的求解结果
    y=x1(2);                      % y的求解结果
    z=x1(3);                      % z的求解结果
    eval(f);
    F1(i)=x;
    F2(i)=y;
    delta(i)=z;
    
    if F1(i)<0 | F1(i)==0 | F2(i)<0 | F2==0
        disp('F出现为0的情况')
    end
end