# Numerical Solvers

### Matthew Kafker, Jeremy Welsh-Kavan

<hr>

The goal of this project is to build a program for numerically solving various ODEs and PDEs. By validating various numerical methods on DEs with known analytical solutions, we hope to be able to implement these methods for a much larger set DEs.

We start with a naive approach and propagate our solution in time using information from the spatial derivative. This is analogous to the the initial value problem (IVP) in which we hope to solve for $y(t)$ given

$
y'(t) = f(t, y(t)) ~~~~~ y(t_0) = y_0
$

Instead, we are solving for $u(x,t)$ given

$
\partial_t u(x,t) = \partial_{xx} u(x,t) \\
u(x,0) = f(x) \\
u(0,t) = u(L,t) = 0 \\
$

for $x \in [0,L]$ and $t>0$.


[//]: <> (For instructions on how to use LaTeX, see the LaTeX folder, which also contains a few other TeXShop files.)

<hr>

Most of what we implement here will be informed by the following resources:

Asmar, Nakhlé H., and Nakhlé H. Asmar. Partial Differential Equations with Fourier Series and Boundary Value Problems. 2nd ed. Upper Saddle River, N. j: Pearson Prentice Hall, 2005. Print.

Press, William H. Numerical Recipes : the Art of Scientific Computing. Cambridge ;: Cambridge University Press, 1989. Print.

Kutz, Jose Nathan. Data-Driven Modeling & Scientific Computation : Methods for Complex Systems & Big Data. First edition. Oxford: Oxford University Press, 2013. Print.