function test(p,n)
x=1:0.01:2;
f1=(x.^4)./5+1./(5.*x);
s=0;
for i=1:n
    s=s+p(i,2)*sigmoid(-p(i,1).*x+p(i,3));
end
f2=2/5+(x-1).*s;
plot(x,f1);
hold on
plot(x,f2);
end