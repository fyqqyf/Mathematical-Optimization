function [p,val,k]=bfgs(fun,gfun,m,n,varargin)
%����: ��BFGS�㷨�����Լ������: min f(p)
%����: p0�ǳ�ʼ��, fun, gfun�ֱ���Ŀ�꺯�������ݶ�;
% varargin������Ŀɱ��������, �򵥵���bfgsʱ���Ժ�����,
% ������������ѭ�����øó���ʱ��������Ҫ������
%mΪѵ������Ŀ��nΪ��Ԫ��Ŀ
%���: p, val�ֱ��ǽ������ŵ������ֵ, k�ǵ�������.
maxk=10000; %��������������
x=linspace(1,2,m);%����ѵ����
p0=rand(n,3);%������ʼֵ����randΪ��ʼֵ)
rho=0.5;sigma=0.5; epsilon=1e-1;
k=0;
Bk=ones(n,n,3);
Bk(:,:,1)=eye(n);%Bk=feval(��Hess��,x0);
Bk(:,:,2)=eye(n);
Bk(:,:,3)=eye(n);
while(k<maxk)
    gk=feval(gfun,x,p0,n,m,varargin{:}); %�����ݶ�
    if(norm(gk)<epsilon), break; end %������ֹ׼��
    %gk(:,1)=gk(:,1)/norm(gk(:,1));
    %gk(:,2)=gk(:,2)/norm(gk(:,2));
    %gk(:,3)=gk(:,3)/norm(gk(:,3));
    a=-Bk(:,:,1)\gk;
    b=-Bk(:,:,2)\gk;
    c=-Bk(:,:,3)\gk;
    dk(:,1)=a(:,1); %�ⷽ����, ������������(ע�������!)
    dk(:,2)=b(:,2);
    dk(:,3)=c(:,3);
    %dk=dk./norm(dk(1,:));
    m0=0; mk=0;
    while(m0<20) % ��Armijo�����󲽳�
        newf=feval(fun,x,p0+rho^m0*dk,n,m,varargin{:});
        oldf=feval(fun,x,p0,n,m,varargin{:});
        newg=feval(gfun,x,p0+rho^m0*dk,n,m);
        if(newf<oldf+sigma*rho^m0*gk'*dk)
            mk=m0; break;
        end
        m0=m0+1;
    end
    %BFGSУ��
    p=p0+rho^mk*dk;
    sk=p-p0; yk=feval(gfun,x,p,n,m,varargin{:})-gk;
    %yk(:,1)=yk(:,1)/norm(yk(:,1));
    %yk(:,2)=yk(:,2)/norm(yk(:,2));
    %yk(:,3)=yk(:,3)/norm(yk(:,3));
    if(yk'*sk>0)
        Bk(:,:,1)=Bk(:,:,1)-(Bk(:,:,1)*sk(:,1)*sk(:,1)'*Bk(:,:,1))/(sk(:,1)'*Bk(:,:,1)*sk(:,1))+(yk(:,1)*yk(:,1)')/(yk(:,1)'*sk(:,1));
        Bk(:,:,2)=Bk(:,:,2)-(Bk(:,:,2)*sk(:,2)*sk(:,2)'*Bk(:,:,2))/(sk(:,2)'*Bk(:,:,2)*sk(:,2))+(yk(:,2)*yk(:,2)')/(yk(:,2)'*sk(:,2));
        Bk(:,:,3)=Bk(:,:,3)-(Bk(:,:,3)*sk(:,3)*sk(:,3)'*Bk(:,:,3))/(sk(:,3)'*Bk(:,:,3)*sk(:,3))+(yk(:,3)*yk(:,3)')/(yk(:,3)'*sk(:,3));
    end
    k=k+1; p0=p;
end
val=feval(fun,x,p0,n,m,varargin{:});