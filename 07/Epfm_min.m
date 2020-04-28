function [ x,y ] = Epfm_min( fx,gx,hx,xx0,s,c,a)
%fx是目标函数
%gx是不等式约束方程组(且g>=0)
%xx0是初始点
%hx是等式约束方程组(且h=0)
%s是精确度(s>0)
%c是放大系数(c>1)
%a是罚因子(默认为1)
syms x1 x2 x3
xx1=xx0;
v=xx0;
a1=a;
Pxk=1;
G=-feval(gx,v);%用于判别max{0,-g(x)}
while Pxk>s
if(G<0)
    Px=a1*hx*hx;
else
    Px=a1*hx*hx+a1*gx*gx;
end
Fx=fx+a1*Px;%将约束问题化为了一个无约束的问题
% 接下来解min F（x）
dFx1=diff(Fx,x1);%分别对x1，x2求偏导数
dFx2=diff(Fx,x2);
dFx3=diff(Fx,x3);
[k,b]=solve(dFx1,dFx2,dFx3,'x1','x2','x3');%求出
xx2=xx1+[k,b];
Pxk=a1*subs(Px,v,xx2);
xx1=xx2;%相当于置k=k+1
a1=c*a1;%罚因子放大
G=-subs(gx,v,xx1);%用于判别max{0,-g(x)}

end
x=xx1;
y=a1/c;