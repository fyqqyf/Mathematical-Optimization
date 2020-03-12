function He=Hess1(x)
n=length(x);
He=zeros(n,n);
He=[8-2*x(2),-2*x(1);
    -2*x(1),2];