function [x,val,k]=dfp(fun,gfun,x0)
%功能: 用DFP算法求解无约束问题: min f(x)
%输入: x0是初始点, fun, gfun分别是目标函数及其梯度
%输出: x, val分别是近似最优点和最优值, k是迭代次数.
%m为训练集数目，n为神经元数目
maxk=1e5; %给出最大迭代次数
rho=0.55;sigma=0.4; epsilon=1e-5;
k=0; n=length(x0);
x=linspace(1,2,m);%构建训练集
p0=ones(n,3);%构建初始值（以rand为初始值)
Hk=ones(n,n,3);
Hk(:,:,1)=5.*eye(n);%Bk=feval(’Hess’,x0);
Hk(:,:,2)=5.*eye(n);
Hk(:,:,3)=5.*eye(n);
dk=zeros(n,3);
while(k<maxk)
gk=feval(gfun,x,p0,n,m,varargin{:}); %计算梯度
if(norm(gk)<epsilon), break; end %检验终止准则
    a=-Hk(:,:,1)\gk;
    b=-Hk(:,:,2)\gk;
    c=-Hk(:,:,3)\gk;
    dk(:,1)=a(:,1); %解方程组, 计算搜索方向(注意是左乘!)
    dk(:,2)=b(:,2);
    dk(:,3)=c(:,3);
    dk=dk/norm(dk(1,:));
dk=-Hk*gk; %解方程组, 计算搜索方向
m=0; mk=0;
while(m<20) % 用Armijo搜索求步长
if(feval(fun,x0+rho^m*dk)<feval(fun,x0)+sigma*rho^m*gk'*dk)
mk=m; break;
end
m=m+1;
end
%DFP校正
x=x0+rho^mk*dk;
sk=p-p0; yk=feval(gfun,x,p,n,m,varargin{:})-gk;
 yk(1,:)=yk(1,:)/norm(yk(1,:));
    yk(2,:)=yk(2,:)/norm(yk(2,:));
    yk(3,:)=yk(3,:)/norm(yk(3,:));
if(sk'*yk>0)
 Hk(:,:,1)=Hk(:,:,1)-(Hk(:,:,1)*sk(:,1)*sk(:,1)'*Hk(:,:,1))/(sk(:,1)'*Hk(:,:,1)*sk(:,1))+(yk(:,1)*yk(:,1)')/(yk(:,1)'*sk(:,1));
        Hk(:,:,2)=Hk(:,:,2)-(Hk(:,:,2)*sk(:,2)*sk(:,2)'*Hk(:,:,2))/(sk(:,2)'*Hk(:,:,2)*sk(:,2))+(yk(:,2)*yk(:,2)')/(yk(:,2)'*sk(:,2));
        Hk(:,:,3)=Hk(:,:,3)-(Hk(:,:,3)*sk(:,3)*sk(:,3)'*Hk(:,:,3))/(sk(:,3)'*Hk(:,:,3)*sk(:,3))+(yk(:,3)*yk(:,3)')/(yk(:,3)'*sk(:,3));
end
k=k+1; x0=x;
end
val=feval(fun,x,p0,n,m,varargin{:});