function [x,val,k]=dampnm(fun,gfun,Hess,x0)
%功能: 用阻尼牛顿法求解无约束问题: min f(x)
%输入: x0是初始点, fun, gfun, Hess 分别是求
% 目标函数值,梯度,Hesse 阵的函数
%输出: x, val分别是近似最优点和最优值, k是迭代次数.
maxk=100; %给出最大迭代次数
rho=0.55;beta=0.5;
k=0; epsilon=1e-3;
figure(1)
clf
hold on
fsurf(@(x,y) 4*(x).^2+(y).^2-8*(x)-4*(y),'EdgeColor','c',...
    'ShowContours','on',...
    'FaceAlpha',0.5);
plot3(x0(1),x0(2),feval(fun,x0),'.','LineWidth',1,...
    'MarkerEdgeColor','y',...
    'MarkerFaceColor','y',...
    'MarkerSize',20)
while(k<maxk)
    gk=feval(gfun,x0); %计算梯度
    Gk=feval(Hess,x0); %计算Hesse阵
    dk=-Gk\gk; %解方程组Gk*dk=-gk, 计算搜索方向
    if(norm(gk)<epsilon), break; end %检验终止准则
    m=0; mk=0;
    while(m<20) % 用Armijo搜索求步长
        if(feval(fun,x0+rho^m*dk)<feval(fun,x0)+beta*rho^m*gk'*dk)
            mk=m; break;
        end
        m=m+1;
    end
    x1=x0+rho^mk*dk;
    plot3(x1(1),x1(2),feval(fun,x1),'.','LineWidth',1,...
    'MarkerEdgeColor','y',...
    'MarkerFaceColor','y',...
    'MarkerSize',20)
    p=plot3([x0(1);x1(1)],[x0(2),x1(2)],[feval(fun,x0),feval(fun,x1)],'b');
    x0=x1;
    k=k+1;
    
end
x=x0;
val=feval(fun,x);
plot3(x(1),x(2),val,'.','LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','k',...
    'MarkerSize',20)