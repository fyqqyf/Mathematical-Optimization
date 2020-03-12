function Armijio(xk,dk)
rho=0.6;%题目中rho=0.6
beta=0.5;%取beta为0.5
m=0; mmax=70;%实验证明只需69次即可收敛
if(fun(xk+dk)<=fun(xk)+rho*gfun(xk)'*dk)
    alpha=1;
else
    while (m<=mmax)
        if(fun(xk+beta*rho^m*dk)<=fun(xk)+beta*rho^(m+1)*gfun(xk)'*dk)
            mk=m;
            break;
        end
        m=m+1;
    end
    alpha=beta*rho^mk;
end

alpha
beta
xk1=xk+alpha*dk
fk=fun(xk)
fk1=fun(xk1)