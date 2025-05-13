## Lid Driven Cavity(顶盖方腔流)

![Schematic](D:\Desktop\latticeboltzmann-python\Lid-Driven Cavity.png)

### 模拟条件



1. 正方形边长 $L=0.2m$

2. 运动粘度 $\upsilon =1.2\times 10^{-3} m^{-2}s^{-1}$

3. 顶盖运动速度 $U=6 m/s$

我那么可以得到 $Re = \frac{L U \rho }{\mu} = \frac{L U}{ \upsilon} = 0.2*6/1.2*10^{-3}=1000$

我们的LBM 计算主要是上述几个变量



### 参数映射

如何从连续模型to 离散模型，个人认为是LBM中，最让人迷惑的一点。

引用书中的话“在LBM中，只要满足雷诺数$Re=1000$就可以任意选择$U$和$\upsilon$的值，为了减少LBM的压缩效应。”

取$U=0.1$,$\upsilon=0.01$ ,$\delta t=1,dx=dy=1$,能够算出 格子数$N=100$



### LBM 基础知识

LBM，格子玻尔兹曼方法，是一种解算玻尔兹曼方程的方法。 
$$
\frac{\partial f}{\partial t}+\mathbf{c} \cdot \nabla f = \Omega
$$
带源项的对流方程，方程左侧表示对流迁移，方程右侧表示碰撞。

我们略去复杂的原理部分，直接给出解算LBM程序框架
$$
f_i(\mathbf{r}+c_i\delta_t,t+\delta_t) - f_i(\mathbf{r},t) = \Omega_i(\mathbf{r},t)
$$
上式中$i = 1...M$ 为离散的格点形式，$c_i=(0,0)...(1,1)$是中心格点指向其可迁移格点的速度大小。如图所示是D2Q9离散格式

![D2Q9](D:\Desktop\latticeboltzmann-python\D2Q9.png)

对碰撞项使用BGK 或 BGKW  近似
$$
\Omega = \omega(f^{\rm{eq}}-f),\omega = 1/\tau
$$

$$
f_i(\mathbf{r}+\mathbf{c}_i\delta t,t+\delta t) = f_i(\mathbf{r},t)+\omega \delta t [f_i^{\rm{eq}}(\mathbf{r},t)-f_i(\mathbf{r},t) ]
$$

使用C-E展开，可以建立运动粘度$\upsilon$和松弛频率$\omega$的关系:
$$
\upsilon = c_s^2\delta t(\frac{1}{\omega}-\frac{1}{2}) \Rightarrow \omega =\frac{1}{3\upsilon+0.5}
$$
$c_s$一般被称为声速，本文中：$c_s=\frac{c}{\sqrt 3}$。

给出一些宏观量 和 LBM中的一些量的联系：
$$
\sum_{i=0}^{8} f_i(\mathbf{r}) = \rho(\mathbf{r})
$$

$$
\sum_{i=0}^{8} \mathbf{c_i}f_i(\mathbf{r}) = \rho(\mathbf{r})\mathbf{u}(\mathbf{r})
$$



一般情况下，在程序中把主方程分成两步（碰撞和迁移）。其中迁移需要结合边界条件，而碰撞则需要结合平衡态速度分布函数。

1. 碰撞
   $$
   f = (1-\omega) f+\omega f^{\mathrm{eq}}
   $$
   

   碰撞中最重要的是平衡态函数$f^{\mathrm{eq}}$,我们不加证明的给出常用的平衡态函数。
   $$
   f^{\mathrm{eq}}_i = \xi_i \rho \left[1+\frac{\mathbf{c}_i\cdot \mathbf{u}}  {c_s^2} +\frac{(\mathbf{c}_i\cdot\mathbf{u})^2}{2c_s^4} -\frac{\mathbf{u}^2}{2 c_s^2}\right]
   $$
   式中$\xi_i$是权函数。

2. 迁移

   对于非边界的部分有：
   $$
   f_i(\mathbf{r}+\mathbf{c}_i\delta t,t+\delta t) = f_i(\mathbf{r},t)
   $$
   可以在程序中非常简单的实现。

   边界处的处理，我们放在下一部分



### 边界条件与初始条件

- 反弹边界

  ![image-20250513151225560](C:\Users\Administrator\AppData\Roaming\Typora\typora-user-images\image-20250513151225560.png)

以本图为例，处于边界处的$f_3,f_6,f_7$都可以通过迁移获得，而$f_1,f_5,f_8$ 则不可以（因为墙里没有格点）。很自然的想法，流向墙里并不符合物理、
$$
f_5(0,t) = f_7(0,t),\ f_1(0,t)=f_3(0,t),\ f_8(0,t)=f_6(0,t)
$$


- 速度边界

  在本小节中$u,v$分别是x方向和y方向上的速度。我们想要得到$f_2,f_5,f_6,\rho$ ，我们已知如下三个方程。
  $$
  \rho = \sum_{i=0}^{8} f_i
  $$

  $$
  \rho u = f_1+f_5+f_8 - f_6-f_3-f_7
  $$

  $$
  \rho v = f_5+f_2+f_6 - f_7-f_4-f_8
  $$

   因此需额外引入一个方程
  $$
  f_2 -f_2^{eq} = f_4 -f_4^{eq}
  $$
  这样即可求解出4个变量
  $$
  f_4 = f_2 - \frac{2}{3}\rho v
  $$

  $$
  f_7 = f_5+\frac{1}{2}(f_1-f_3) - \frac{1}{6}\rho v-\frac{1}{2}\rho u
  $$

  $$
  f_8= f_6 +\frac{1}{2}(f_3-f_1) - \frac{1}{6}\rho v+\frac{1}{2}\rho u 
  $$

  **不过书中所写的代码与 书中给的边界条件不同**, 十分令人疑惑

  将能够产生正确速度分布的代码注记于下：

  ````python
      rhon = F[ 1:-1,-1,0] + F[ 1:-1 ,-1 ,1] + F[1:-1, -1,3]  +2(F[1:-1,-1,2]+F[1:-1,-1,6]+F[1:-1,-1,5])
  
      F[1:-1,-1,4] = F[1:-1,-1,2]
      F[1:-1,-1,8] = F[1:-1,-1,6] + rhon*uo/6
      F[1:-1,-1,7] = F[1:-1,-1,5] - rhon*uo/6
  ````

  

猜测是在边界处应满足$f_3 -f_1 = \frac{2}{3}\rho u $ ，因为这与 $f_2 -f_2^{eq} = f_4 -f_4^{eq}$ 可以得到相同的结果。



