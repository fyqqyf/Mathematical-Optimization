function He=Hess(x)
n=length(x);
He=zeros(n,n);
He=[8,0;
    0,2];