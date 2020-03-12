function [s,phis,k,G,E]=golds(phi,a,b,delta,epsilon)
%����: phi��Ŀ�꺯��, a, b ����������������˵�
% delta, epsilon�ֱ����Ա����ͺ���ֵ���������
% ���: s, phis�ֱ��ǽ��Ƽ�С��ͼ�Сֵ, G��n x 4����,
% ���k�зֱ���a,p,q,b�ĵ�k�ε���ֵ[ak,pk,qk,bk],
% E=[ds,dphi], �ֱ���s��phis�������.
x=0:0.0001:1;
y=feval(phi,x);
plot(x,y);
hold on;
t=(sqrt(5)-1)/2;%������
h=b-a;%��ʼ���䳤��
phia=feval(phi,a); phib=feval(phi,b);%�����ʱ���˵㺯��ֵ

p=a+(1-t)*h; q=a+t*h;
phip=feval(phi,p); phiq=feval(phi,q);

k=1;
G(k,:)=[a, p, q, b];

while(abs(phib-phia)>epsilon)||(h>delta)
    if(phip<phiq)
        b=q; phib=phiq; q=p; phiq=phip;
        h=b-a; p=a+(1-t)*h; phip=feval(phi,p);
    else
        a=p; phia=phip; p=q; phip=phiq;
        h=b-a; q=a+t*h; phiq=feval(phi,q);
    end
    k=k+1; 
    G(k,:)=[a, p, q, b];
    %text(a,phia,'o');
    %text(b,phib,'o');
    plot([a;b],[phia,phib]);
    hold on;
end
ds=abs(b-a); dphi=abs(phib-phia);
if(phip<=phiq)
    s=p; phis=phip;
else
    s=q; phis=phiq;
end
text(s,phis,'*','color','g');
E=[ds,dphi];