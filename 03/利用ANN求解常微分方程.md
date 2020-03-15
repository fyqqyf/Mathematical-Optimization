# 利用ANN求解常微分方程

在《最优化方法》上学习拟牛顿法（BFGS，DFP）后,我们可以尝试使用这两种优化方法进行ANN（人工神经网络）的搭建并且使用它来求解常微分方程。

<img src="C:\Users\50672\AppData\Roaming\Typora\typora-user-images\image-20200315173954019.png" alt="image-20200315173954019" style="zoom:67%;" />

前馈神经网络由输入层，隐藏层，输出层构成。

考虑一节常微分方程

​                                                               $$\frac{dy(x)}{dx}=f(x,y)$$

初始条件为$y(0)=y_0$.利用神经网络我们可以得到一个实验解$y_t(x,p)$,用来近似方程的精确解，其中$p$即为我们需要训练的参数

$y_t(x,p)=y_0+xN(x,p)$

$N(x,p)=\sum_{j=1}^{n}v_j \varphi (z_j) z_j=\omega_j x-\theta_j$

$\varphi(z_j)$为sigmoid函数

![image-20200315175102208](C:\Users\50672\AppData\Roaming\Typora\typora-user-images\image-20200315175102208.png)

使用最小二乘作为loss函数

分别使用BFGS和DFP进行优化

(使用$MATLAB$进行编程，遇到了如归一化，奇异精度等问题，还有$fun$和$gfun$函数的编写，都有很大的改进空间）









