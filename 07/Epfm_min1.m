function [ x,y ] = Epfm_min1( fx,gx,hx,xx0,s,c,a)
%fx是目标函数
%gx是不等式约束方程组(且g>=0)
%xx0是初始点
%hx是等式约束方程组(且h=0)
%s是精确度(s>0)
%c是放大系数(c>1)
%a是罚因子(默认为1)
xx1=xx0;
v=xx0;
a1=a;
Pxk=1;%假设Px等于1，以免不必要错误
while Pxk>s
    g=feval(gx,v);
    h=feval(hx,v);
    f=feval(fx,v);
    G=max(-g,0);%用于判别max{0,-g(x)}
    if(G<0)
        Px=a1*h*h;
    else
        Px=a1*h*h+a1*g*g;
    end
    Fx=f+a1*Px;%将约束问题化为了一个无约束的问题
    % 接下来解min F（x）
    [x,ival,ik]=bfgs('Fx','gFx',x0,fun,hf,gf,dfun,dhf,dgf,s,c,a);
    %[k,b,e]=solve(dFx1,dFx2,dFx3,'x1','x2','x3');%求出
    xx1=x;
    Pxk=a1*feval(hx,v);
    xx1=xx2;%相当于置k=k+1
    a1=c*a1;%罚因子放大
end
x=xx1;
y=a1/c;