function [x,val,k]=dampnm(fun,gfun,Hess,x0)
%����: ������ţ�ٷ������Լ������: min f(x)
%����: x0�ǳ�ʼ��, fun, gfun, Hess �ֱ�����
% Ŀ�꺯��ֵ,�ݶ�,Hesse ��ĺ���
%���: x, val�ֱ��ǽ������ŵ������ֵ, k�ǵ�������.
maxk=100; %��������������
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
    gk=feval(gfun,x0); %�����ݶ�
    Gk=feval(Hess,x0); %����Hesse��
    dk=-Gk\gk; %�ⷽ����Gk*dk=-gk, ������������
    if(norm(gk)<epsilon), break; end %������ֹ׼��
    m=0; mk=0;
    while(m<20) % ��Armijo�����󲽳�
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