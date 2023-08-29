function Pasokh=bisection_method(y,a,b,eps)
syms x
X=(a+b)/2; % X赋初值
fx=subs(y,X); % 将X带入求出y,给fx赋初值
E=abs(fx); %求绝对值，给E赋初值
while E>eps %利用中心差值法求解，eps为误差
    fa=subs(y,a);
    fb=subs(y,b);
    fx=subs(y,X);
    if fa*fx<0
        b=X;
    else
        a=X;
    end
    X=(a+b)/2;
    fx=subs(y,X);
    E=abs(fx);
end
Pasokh=X;