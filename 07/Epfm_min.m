function [ x,y ] = Epfm_min( fx,gx,hx,xx0,s,c,a)
%fx��Ŀ�꺯��
%gx�ǲ���ʽԼ��������(��g>=0)
%xx0�ǳ�ʼ��
%hx�ǵ�ʽԼ��������(��h=0)
%s�Ǿ�ȷ��(s>0)
%c�ǷŴ�ϵ��(c>1)
%a�Ƿ�����(Ĭ��Ϊ1)
syms x1 x2 x3
xx1=xx0;
v=xx0;
a1=a;
Pxk=1;
G=-feval(gx,v);%�����б�max{0,-g(x)}
while Pxk>s
if(G<0)
    Px=a1*hx*hx;
else
    Px=a1*hx*hx+a1*gx*gx;
end
Fx=fx+a1*Px;%��Լ�����⻯Ϊ��һ����Լ��������
% ��������min F��x��
dFx1=diff(Fx,x1);%�ֱ��x1��x2��ƫ����
dFx2=diff(Fx,x2);
dFx3=diff(Fx,x3);
[k,b]=solve(dFx1,dFx2,dFx3,'x1','x2','x3');%���
xx2=xx1+[k,b];
Pxk=a1*subs(Px,v,xx2);
xx1=xx2;%�൱����k=k+1
a1=c*a1;%�����ӷŴ�
G=-subs(gx,v,xx1);%�����б�max{0,-g(x)}

end
x=xx1;
y=a1/c;