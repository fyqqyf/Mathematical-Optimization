function f=fun(x)
f=0;
for i=1:1e3
    f=f+exp(x(i,1))-x(i,1);
end
