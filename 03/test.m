function test(p,n)
syms x;
f1=(x.^4)./5+1./(5.*x);
s=0;
for i=1:n
    s=s+p(i,2).*sigmoid(-p(i,1).*x+p(i,3));
end
f2=2/5+(x-1).*s;

f3=diff(f1,x);
f4=diff(f2,x);
y=1:0.1:2;
figure(1)
subplot(1,2,1)
plot(y,eval(subs(f1,x,y)));
hold on
plot(y,eval(subs(f2,x,y)));
subplot(1,2,2)
plot(y,eval(subs(f3,x,y)));
hold on 
plot(y,eval(subs(f4,x,y)));
end