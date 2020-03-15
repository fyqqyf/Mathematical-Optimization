function [p,val,k]=bfgs(fun,gfun,m,n,varargin)
%����: ��BFGS�㷨�����Լ������: min f(p)
%����: p0�ǳ�ʼ��, fun, gfun�ֱ���Ŀ�꺯�������ݶ�;
% varargin������Ŀɱ��������, �򵥵���bfgsʱ���Ժ�����,
% ������������ѭ�����øó���ʱ��������Ҫ������
%mΪѵ������Ŀ��nΪ��Ԫ��Ŀ
%���: p, val�ֱ��ǽ������ŵ������ֵ, k�ǵ�������.
maxk=1000; %��������������
x=linspace(1,2,m);%����ѵ����
p0=ones(n,3);%������ʼֵ����randΪ��ʼֵ)
rho=0.5;sigma=0.3; epsilon=1;sita=0.5;
k=0;
%n=length(p0);
Bk=ones(n,n,3);
Bk(:,:,1)=5.*eye(n);%Bk=feval(��Hess��,x0);
Bk(:,:,2)=5.*eye(n);
Bk(:,:,3)=5.*eye(n);
dk=zeros(n,3);
while(k<maxk)
    gk=feval(gfun,x,p0,n,m,varargin{:}); %�����ݶ�
    gk=gk/norm(gk(1,:));
    if(norm(gk)<epsilon), break; end %������ֹ׼��
    a=-Bk(:,:,1)\gk;
    b=-Bk(:,:,2)\gk;
    c=-Bk(:,:,3)\gk;
    dk(:,1)=a(:,1); %�ⷽ����, ������������(ע�������!)
    dk(:,2)=b(:,2);
    dk(:,3)=c(:,3);
    dk=dk/norm(dk(1,:));
    m0=0; mk=0;
    while(m0<20) % ��Armijo�����󲽳�
        newf=feval(fun,x,p0+rho^m0*dk,n,m,varargin{:});
        %newg=feval(gfun,x,p0+rho^m0*dk,n,m,varargin{:});
        oldf=feval(fun,x,p0,n,m,varargin{:});
        if(newf<oldf+sigma*rho^m0*gk'*dk)
            mk=m0; break;
        end
        m0=m0+1;
    end
    %BFGSУ��
    p=p0+rho^mk*dk;
    sk=p-p0; yk=feval(gfun,x,p,n,m,varargin{:})-gk;
    yk(1,:)=yk(1,:)/norm(yk(1,:));
    yk(2,:)=yk(2,:)/norm(yk(2,:));
    yk(3,:)=yk(3,:)/norm(yk(3,:));
    if(norm(yk'*sk)>0)
        Bk(:,:,1)=Bk(:,:,1)-(Bk(:,:,1)*sk(:,1)*sk(:,1)'*Bk(:,:,1))/(sk(:,1)'*Bk(:,:,1)*sk(:,1))+(yk(:,1)*yk(:,1)')/(yk(:,1)'*sk(:,1));
        Bk(:,:,2)=Bk(:,:,2)-(Bk(:,:,2)*sk(:,2)*sk(:,2)'*Bk(:,:,2))/(sk(:,2)'*Bk(:,:,2)*sk(:,2))+(yk(:,2)*yk(:,2)')/(yk(:,2)'*sk(:,2));
        Bk(:,:,3)=Bk(:,:,3)-(Bk(:,:,3)*sk(:,3)*sk(:,3)'*Bk(:,:,3))/(sk(:,3)'*Bk(:,:,3)*sk(:,3))+(yk(:,3)*yk(:,3)')/(yk(:,3)'*sk(:,3));
    end
    k=k+1; p0=p;
end
val=feval(fun,x,p0,n,m,varargin{:});