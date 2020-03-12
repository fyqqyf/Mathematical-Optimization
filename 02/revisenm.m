function [x,val,k]=revisenm(fun,gfun,Hess,x0)
% ����: ������ţ�ٷ������Լ������: min f(x)
%����: x0�ǳ�ʼ��, fun, gfun, Hess �ֱ�����
% Ŀ�꺯��ֵ,�ݶ�,Hesse ��ĺ���
%���: x, val�ֱ��ǽ������ŵ������ֵ, k�ǵ�������.
n=length(x0); maxk=150;
rho=0.55;sigma=0.5; tau=0.0;
k=0; epsilon=1e-3;
figure(1)
clf
hold on
fsurf(@(x,y) 4*x^2+y^2-x^2*y,'EdgeColor','none',...
    'ShowContours','on',...
    'FaceAlpha',0.5);
plot3(x0(1),x0(2),feval(fun,x0),'.','LineWidth',1,...
    'MarkerEdgeColor','y',...
    'MarkerFaceColor','y',...
    'MarkerSize',20)
while(k<maxk)
    gk=feval(gfun,x0); % �����ݶ�
    muk=norm(gk)^(1+tau);
    Gk=feval(Hess,x0); % ����Hesse��
    Ak=Gk+muk*eye(n);
    dk=-Ak\gk; %�ⷽ����Gk*dk=-gk, ������������
    if(norm(gk)<epsilon), break; end %������ֹ׼��
    m=0; mk=0;
    while(m<20) %��Armijo�����󲽳�
        if(feval(fun,x0+rho^m*dk)<feval(fun,x0)+sigma*rho^m*gk'*dk)
            mk=m; break;
        end
        m=m+1;
    end
    x1=x0+rho^mk*dk;
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