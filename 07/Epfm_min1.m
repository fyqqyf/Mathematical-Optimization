function [ x,y ] = Epfm_min1( fx,gx,hx,xx0,s,c,a)
%fx��Ŀ�꺯��
%gx�ǲ���ʽԼ��������(��g>=0)
%xx0�ǳ�ʼ��
%hx�ǵ�ʽԼ��������(��h=0)
%s�Ǿ�ȷ��(s>0)
%c�ǷŴ�ϵ��(c>1)
%a�Ƿ�����(Ĭ��Ϊ1)
xx1=xx0;
v=xx0;
a1=a;
Pxk=1;%����Px����1�����ⲻ��Ҫ����
while Pxk>s
    g=feval(gx,v);
    h=feval(hx,v);
    f=feval(fx,v);
    G=max(-g,0);%�����б�max{0,-g(x)}
    if(G<0)
        Px=a1*h*h;
    else
        Px=a1*h*h+a1*g*g;
    end
    Fx=f+a1*Px;%��Լ�����⻯Ϊ��һ����Լ��������
    % ��������min F��x��
    [x,ival,ik]=bfgs('Fx','gFx',x0,fun,hf,gf,dfun,dhf,dgf,s,c,a);
    %[k,b,e]=solve(dFx1,dFx2,dFx3,'x1','x2','x3');%���
    xx1=x;
    Pxk=a1*feval(hx,v);
    xx1=xx2;%�൱����k=k+1
    a1=c*a1;%�����ӷŴ�
end
x=xx1;
y=a1/c;