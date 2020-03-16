function [p,val,k]=dfp(fun,gfun,m,n)
%����: ��DFP�㷨�����Լ������: min f(x)
%����: x0�ǳ�ʼ��, fun, gfun�ֱ���Ŀ�꺯�������ݶ�
%���: x, val�ֱ��ǽ������ŵ������ֵ, k�ǵ�������.
%mΪѵ������Ŀ��nΪ��Ԫ��Ŀ
maxk=1e5; %��������������
rho=0.5;sigma=0.5; epsilon=1e-1;
k=0;
x=linspace(1,2,m);%����ѵ����
p0=rand(n,3);%������ʼֵ����randΪ��ʼֵ)
Hk=ones(n,n,3);
Hk(:,:,1)=eye(n);%Bk=feval(��Hess��,x0);
Hk(:,:,2)=eye(n);
Hk(:,:,3)=eye(n);
dk=zeros(n,3);
while(k<maxk)
    gk=feval(gfun,x,p0,n,m); %�����ݶ�
    if(norm(gk)<epsilon), break; end %������ֹ׼��
    a=-Hk(:,:,1)*gk;
    b=-Hk(:,:,2)*gk;
    c=-Hk(:,:,3)*gk;
    dk(:,1)=a(:,1); %�ⷽ����, ������������(ע�������!)
    dk(:,2)=b(:,2);
    dk(:,3)=c(:,3);
    %dk=dk/norm(dk(1,:));
    %dk=-Hk*gk; %�ⷽ����, ������������
    m0=0; mk=0;
    while(m0<20) % ��Armijo�����󲽳�
        if(feval(fun,x,p0+rho^m0*dk,n,m)<feval(fun,x,p0,n,m)+sigma*rho^m0*gk'*dk)
            mk=m0; break;
        end
        m0=m0+1;
    end
    %DFPУ��
    p=p0+rho^mk*dk;
    sk=p-p0; yk=feval(gfun,x,p,n,m)-gk;
    %yk(1,:)=yk(1,:)/norm(yk(1,:));
    %yk(2,:)=yk(2,:)/norm(yk(2,:));
    %yk(3,:)=yk(3,:)/norm(yk(3,:));
    if(sk'*yk>0)
        Hk(:,:,1)=Hk(:,:,1)-(Hk(:,:,1)*sk(:,1)*sk(:,1)'*Hk(:,:,1))/(sk(:,1)'*Hk(:,:,1)*sk(:,1))+(yk(:,1)*yk(:,1)')/(yk(:,1)'*sk(:,1));
        Hk(:,:,2)=Hk(:,:,2)-(Hk(:,:,2)*sk(:,2)*sk(:,2)'*Hk(:,:,2))/(sk(:,2)'*Hk(:,:,2)*sk(:,2))+(yk(:,2)*yk(:,2)')/(yk(:,2)'*sk(:,2));
        Hk(:,:,3)=Hk(:,:,3)-(Hk(:,:,3)*sk(:,3)*sk(:,3)'*Hk(:,:,3))/(sk(:,3)'*Hk(:,:,3)*sk(:,3))+(yk(:,3)*yk(:,3)')/(yk(:,3)'*sk(:,3));
    end
    k=k+1; p0=p;
end
val=feval(fun,x,p0,n,m);