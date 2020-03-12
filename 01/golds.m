function [s,phis,k,G,E]=golds(phi,a,b,delta,epsilon)
%输入: phi是目标函数, a, b 是搜索区间的两个端点
% delta, epsilon分别是自变量和函数值的容许误差
% 输出: s, phis分别是近似极小点和极小值, G是n x 4矩阵,
% 其第k行分别是a,p,q,b的第k次迭代值[ak,pk,qk,bk],
% E=[ds,dphi], 分别是s和phis的误差限.
x=0:0.0001:1;
y=feval(phi,x);
plot(x,y);
hold on;
t=(sqrt(5)-1)/2;%缩短率
h=b-a;%初始区间长度
phia=feval(phi,a); phib=feval(phi,b);%计算此时两端点函数值

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