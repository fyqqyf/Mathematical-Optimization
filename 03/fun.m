function f=fun(x,p,n,m)
f=0;
for i=1:m
    s1=0;
    s2=0;
    for j=1:n
        x0=-p(j,1).*x(i)+p(j,3);
        s0=sigmoid(x0);
        s1=s1+p(j,2).*s0;
        s2=s2+s0^2.*p(j,1).*p(j,2).*exp(x0);
    end
    f0=(2-1/x(i)).*s1+(x(i)-1).*s2-x(i)^3+(2/5)/x(i);
    f=f+(f0).^2;
end
end



