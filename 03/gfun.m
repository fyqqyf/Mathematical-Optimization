function f=gfun(x,p,n,m)
f=0;


%syms p1;
%syms p2;
%syms p3;
%syms x0;
%f1(p1,p2,p3,x0)=p2*sigmoid(-p1*x0+p3);
%f11(p1,p2,p3,x0)=diff(f1,p1);
%f12(p1,p2,p3,x0)=diff(f1,p2);
%f13(p1,p2,p3,x0)=diff(f1,p3);
%f2(p1,p2,p3,x0)=exp(-p1.*x0+p3).*p1.*p2.*(sigmoid(-p1*x0+p3)^2);
%f21(p1,p2,p3,x0)=diff(f1,p1);
%f22(p1,p2,p3,x0)=diff(f1,p2);
%f23(p1,p2,p3,x0)=diff(f1,p3);
%f11(p1,p2,p3,x0)=(p2*x0*exp(p3 - p1*x0))/(exp(p3 - p1*x0) + 1)^2;
%f12(p1, p2, p3, x0) =1/(exp(p3 - p1*x0) + 1);
%f13(p1, p2, p3, x0) =-(p2*exp(p3 - p1*x0))/(exp(p3 - p1*x0) + 1)^2;
%f21(p1, p2, p3, x0) =(p2*exp(p3 - p1*x0))/(exp(p3 - p1*x0) + 1)^2 - (p1*p2*x0*exp(p3 - p1*x0))/(exp(p3 - p1*x0) + 1)^2 + (2*p1*p2*x0*exp(2*p3 - 2*p1*x0))/(exp(p3 - p1*x0) + 1)^3
%f22(p1, p2, p3, x0) = (p1*exp(p3 - p1*x0))/(exp(p3 - p1*x0) + 1)^2
%f23(p1, p2, p3, x0) =(p1*p2*exp(p3 - p1*x0))/(exp(p3 - p1*x0) + 1)^2 - (2*p1*p2*exp(2*p3 - 2*p1*x0))/(exp(p3 - p1*x0) + 1)^3
%需要时可用 subs eval 求值（但比较慢）
for i=1:m
    s1=0;
    s2=0;
    v1=zeros(n,3);
    v2=zeros(n,3);
    for j=1:n
        x1=-p(j,1).*x(i)+p(j,3);
        s0=sigmoid(x1);
        s1=s1+p(j,2).*s0;
        s2=s2+s0^2.*p(j,1).*p(j,2).*exp(x1);
        v1(j,1)=(p(j,2)*x(i)*exp(p(j,3) - p(j,1)*x(i)))/(exp(p(j,3) - p(j,1)*x(i)) + 1)^2;
        v1(j,2)=1/(exp(p(j,3) - p(j,1)*x(i)) + 1);
        v1(j,3)=-(p(j,2)*exp(p(j,3) - p(j,1)*x(i)))/(exp(p(j,3) - p(j,1)*x(i)) + 1)^2;
        %v1=v1+[f14,f15,f16]';
        v2(j,1)=(p(j,2)*exp(p(j,3) - p(j,1)*x(i)))/(exp(p(j,3) - p(j,1)*x(i)) + 1)^2 - (p(j,1)*p(j,2)*x(i)*exp(p(j,3) - p(j,1)*x(i)))/(exp(p(j,3) - p(j,1)*x(i)) + 1)^2 + (2*p(j,1)*p(j,2)*x(i)*exp(2*p(j,3) - 2*p(j,1)*x(i)))/(exp(p(j,3) - p(j,1)*x(i)) + 1)^3;
        v2(j,2)=(p(j,1)*exp(p(j,3) - p(j,1)*x(i)))/(exp(p(j,3) - p(j,1)*x(i)) + 1)^2;
        v2(j,3)=(p(j,1)*p(j,2)*exp(p(j,3) - p(j,1)*x(i)))/(exp(p(j,3) - p(j,1)*x(i)) + 1)^2 - (2*p(j,1)*p(j,2)*exp(2*p(j,3) - 2*p(j,1)*x(i)))/(exp(p(j,3) - p(j,1)*x(i)) + 1)^3;
        %v2=v2+[f24,f25,f26]';
    end
    v3=(2-1/x(i)).*v1+(x(i)-1).*v2;
    f0=(2-1/x(i)).*s1+(x(i)-1).*s2-x(i)^3+(2/5)/x(i);
    f=f+2.*f0.*v3;
end
%f=f';
end