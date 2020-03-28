function qd=qk(x,d)
gk=gfun(x);gk=gk';
Bk=Hess(x);
qd=gk'*d+0.5*d'*Bk*d;