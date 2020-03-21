function g=gfun(x)
g=zeros(1e3,1);
for i=1:1e3
    g(i,1)=exp(x(i,1))-1;
end