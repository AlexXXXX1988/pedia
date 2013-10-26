核密度估计(kernel density estimation)
===================================

_分类：非参数统计_

概述
---

_统计学_中，核密度估计是一种_非参数统计_方法。核函数估计用于对样本的频率分布进行平滑变换，近似估计总体的_概率密度函数_。该方法由Rosenblatt(1955)和Emanuel Parzen(1962)共同提出,又名Parzen窗(Parzen window)。

定义
---

设`$(x_1,x_2,...,x_n)$`为取自服从概率密度函数_f_的_独立同分布_样本。为通过样本的频率分布估计总体概率密度函数_f_，为每个样本点选取对称的_核函数_`$K(·)$`，则总体的核函数估计可以表示为

`$\hat{f}_h(x)=\frac{1}{n}\sum_{i=1}{n}K_h(x-x_i)=\frac{1}{nh}\sum_{i=1}{n}K(\frac{x-x_i}{h})$`

上式中，`$K(·)$`是一个对称函数，定义域内不要求均为正直，求积分和为1。`$h$`为_带宽_，是一个调整平滑程度的参数，`$h$`越大则核密度估计曲线越光滑。`$h$`作为`$K(·)$`的参数有如下的关系，`$K_h(x)=\frac{1}{h}K(\frac{x}{h})$`。 核函数有一系列对称函数可以选择，如_均匀分布核函数，三角核函数，Epanechnikov核函数，biweight核函数，triweight核函数，正态分布核函数_等。除均匀分布核函数意外其余上述核函数都体现了样本点在观测值附近对总体概率密度函数的影响是往两端逐渐减弱。而核密度估计的基本思想是，样本作为一个随机变量，出现在观测值附近的可能性由`$K(·)$`来描述，核密度估计曲线是所有核函数的加总。核函数的把样本点对概率密度函数的影响从点变成域，从而可以平滑样本的频率分布，渐近估计总体概率密度函数。带宽`$h$`描述一个样本点的影响范围，`$h$`越大，则样本点的影响向两端削减越慢，曲线越平滑。
上面的方法对于多元数据同样适用，这时仅需要使用所谓的multivariate kernel函数即可，通常得到如下的估计量
`$\hat{f_H}(x)=\frac{1}{n}\sum_{i=1}{n}\frac{1}{|\textbf{H}|}K(\textbf{H}^{-1}(\textbf{x}-\textbf{x}_i))$`.

我们一般都会采用简单的带宽阵或K(·),比如`$\textbf{H}=diag\{h_1,...,h_d\}$`

`$\hat{f_H}(x)=\frac{1}{nh_1...h_d}\sum_{i=1}{n}\prod{j=1}{d}K(\textbf{H}^{-1}(\textbf{x}-\textbf{x}_i))$`

带宽的选取
--------
我们必须首先考虑如何评价密度估计量的性质。我们来看均方误差`$MSE(\hat{f_h}(x))$`和所谓的积分均方误差（Mean Integrated Square Error）`$MISE(h)$`。它们的定义是
`$MISE(h)=\intMSE(\hat{f_h}(x))dx=\intvar{\hat{f_h}(x)}$`
其中且`$bias\{\hat{f_h}(x)\}=E\{\hat{f_h}-f(x)\}$`显然，`$MSE(\hat{f_h}(x))$`是对于评价估计`$\hat{f_h}(x)$`好坏的一个逐点准则,而`$MISE(h)$`可看成是在每点`$x$`处对局部均方误差的累积,也就是一个全局的评判准则.

进一步地,我们还需要对核函数和f作一定的假设.设K是对称连续的核密度函数,均值为零,方差`$0<\sigma_K^2<\infty$`.这里f有二阶有界连续导数且`$\[f''(t)^2dt<\infty$`,也就是对其光滑性做一定的假设.假设当`$n\rightarrow\infty$`时`$nh\rightarrow\infty$`,`$h\rightarrow0$`,我们将进一步分析该表达式.要计算中的偏差项,注意到应用变量替换有
`$E\{\hat{f_h}(x)\}=\frac{1}{h}\intK\big(\frac{x-u}{h}\big)f(u)du=\intK(t)(x-ht)dt$`
然后用Taylor级数展开
`$f(x-ht)=f(x)-htf'(x)+\frac{h^2t^2}{2}f''(x)+o(h^2)$`
替换并注意到K关于零点对称可得
`$\{E\{\hat{f_h}(x)\}=f(x)+\frac{h^2\sigma_K^2}{2}f''(x)+o(h^2)$`
因此
`$\big(bias\{\hat{f_h}\}\big)^2=h^4\sigma_K^4[f''(x)]^2/4+o(h^4)$`
且该表达式对x积分可得
`$\int\big(bias\{\hat{f_h}\}\big)^2dx=h^4\sigma_K^4\int[f''(x)]^2dx/4+o(h^4)$`
方差项`$\sigma$`可采用类似的方法
`$var\{\hat{f_h}(x)\}=\frac{1}{n}var\{\frac{1}{h}K\big(\frac{x-x_i}{h}\big)\}
=\frac{1}{nh}\intK^2(t)f(x-ht)dt-\frac{1}{n}[E\{\frac{1}{n}K\big(\frac{x-x_i}{h}\big)\}]
=\frac{1}{nh}\intK^2(t)[f(x)+o(1)]dt-\frac{1}{n}[f(x)+o(1)]^2
=\frac{1}{nh}f(x)\intK^2(x)dx+o(\frac{1}{nh})$`
将其对x积分得
`$\intvar\{\hat{f_h}(x)\}=\frac{1}{nh}\intK^2(x)dx+(\frac{1}{nh})$`
因此`$MISE(h)=AMISE(h)+o(\frac{1}{h})$`
`$AMISE(h)=\frac{K^2(x)dx}{nh}+\frac{h^4\sigma^4_K[f''(x)]^2dx}{4}$`
称作渐进均方积分误差.如果当`$n\rightarrow\infty$`时`$nh\rightarrow\infty$`,`$h\rightarrow0$`,则`$MISE(h)\rightarrow0$`
要关于h最小化`$AMISE(h)$`,我们必须把h设在某个中间值,这样可以避免`$\hat{f_h}$`有过大的偏差（太过光滑）或过大的方差（即过于光滑）.关于`$h$`最小化`$AMISE(h)$`表明最好是精确地平衡`$AMISE(h)$`中偏差项和方差项的阶数.显然最优的带宽是
`$h=\big(\frac{\intK^2(x)dx}{n\sigma^4_K[f''(x)]^2dx}\big)^{1/5}~n^{-1/5}$`
上式说明了核函数估计中带宽`$h$`和样本量`$n$`之间的关系，要想得到最佳的带宽`$h$`就要估计上式中的`$f''(x)^2$`
###代入法（Plug-in）
代入法考虑在最优带宽中使用某适当的估计`$R(f'')$`来代替`$R(\hat{f''})$`.在众多的方法中,最简单且最常用的即是Sheatherand Jones(1991;JRSSB)所提出的`$R(f'')=R(\hat{f''})$`的基于核的估计量为
`$\hat{f}''(x)=\frac{\partial^2}{\partialx^2}{\frac{1}{nh_0}\sum_{i=1}{n}\big(\frac{x-x_i}{h_0}\big)}=\frac{1}{nh_0^3}\sum_{i=1}^{n}L''\big(\frac{x-x_i}{h_0}\big)$`
其中`$h_0$`为带宽,`$L$`为用来估计`$f''$`的核函数（不一定与`$K$`相同）.再对其平方并对x积分后即可得到`$R(\hat{f''})$`
###经验法则（Rule of Thumb）
简便起见,我们定义`$R(g)=g^2(z)dz$`.针对最小化AMISE得到的最优带宽中含有未知量`$R(f'')$`,Silverman提出一种初等的方法,rule of thumb（拇指法则,即根据经验的方法）：把`$f$`用方差和估计方差相匹配的正态密度替换.这就等于用`$R(\phi'')/\hat{\sigma}^5$`估计`$R(f'')$`,其中`$\phi$`为标准正态密度函数.若取K为高斯密度核函数而`$\sigma$`使用样本方差`$\sigma$`,Silverman拇指法则得到
`$h=\big(\frac{4}{3n}\big)^{1/5}\hat{\sigma}$`
但通常推荐考虑半极差（IQR,即上下四分位数的距离）,因为IQR是更加稳健的散度度量.Silverman建议在(7.7)中用`$\tilde{\sigma}=min{\sigma,\frac{IQR}{(\Phi^{-1}(0.75)-\Phi^{-1}(0.25))}\}\approx min{\sigma,\frac{IQR}{1.35}}替换`$\sigma$`,其中`$\Phi$`是标准正态累积分布函数.作为产生近似带宽的一种方法,Silverman的拇指法则是很有效的,这种带宽通常作为更加复杂的plugin方法中非参数估计的带宽

###交互验证法（Cross-Validation）
该方法是一种基于数据（data-driven）方法,其基本思想是考虑最小化积分平方误差`$ISE(h)$`.把积分平方误差重新写成
`$ISE(h)=\int\hat{f^2}(x)dx-2E\{\hat{f}(X)\}+\intf^2(x)dx
=R(\hat{f})-2E\{\hat{f}(X)\}+R(f)$`
该表达式的最后一项`$R(f)$`不依赖于`$f$`,因此其对选择`$h$`没有影响.中间项可以用`$\frac{2}{n}\sum_{i=1}^{n}\hat{f}_{-i}(x_i)$`来估计,其中
`$\hat{f}_{-i}(x_i)=\frac{1}{h(n-1)}\sum_{i=1}^nK\big(\frac{x_i-x_j}{h}\big)$`
表示在`$x_i$`点处核密度估计量用除`$x_i$`外所有数据估计的密度,
因此,通过关于`$h$`最小化
`$CV(h)=R(\hat{f})-\frac{2}{n}\sum_{i=1}^n\hat{f}_{-i}(X_i)$`
通常可得到较好的带宽.事实上,可以证明,这样选出来的`$h$`,不妨记为`$\hat{h}_{CV}=arg\minCV(h)$`,按下面的意义是渐进最优的:当`$n\rightarrow\infty$`时,
`$\frac{ISE(\hat{h}_{CV})}{inf_{h>0}ISE(h)}\rightarrow1$`
