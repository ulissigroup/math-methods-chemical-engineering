#!/usr/bin/env python
# coding: utf-8

# {{ badge }}
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$

# # Second Order Differential Equations

# * Occur in many important applications: fluid mechanics, diffusion, heat transport, statics and circuit theory.\
# In general $F(y'',y',y,x)=0$ where $y'' = \frac{d^2y}{dx^2}$
# * Analytical solutions are possible (usually) only when $2^\circ$ ODE is linear. Exceptions are two special types of $2^\circ$ ODE that can be solved by change of variables to make $1^\circ$ ODE. 

# ## Linear $2^\circ$ ODE:

# \begin{align}
# y'' + p(x) y' + q(x) y = r(x)
# \end{align}
# where $p(x), q(x)$ and $r(x)$ are any functions of $x$.\
# $\rightarrow$ if $y''$ has a coefficient, divide through.
# 
# * Linear $2^\circ$ ODEs come in two flavors:
#   1. $r(x) = 0 \implies$ homogeneous
#   2. $r(x) \neq 0 \implies$ non-homogeneous
# 
# $\rightarrow$ We will begin with homogeneous equations. We will later show that once homogeneous equation has been solved, we can always solve the non-homogeneous one. 
# * Because $2^\circ$ ODEs typically require two integrations, we will have two constants of integration. \
# $\therefore$ We need two initial/boundary conditions to find particular solutions. 
# \begin{align}
# \text{usually} && y(x_0) = k1 &&\text{&} && y'(x_0) = k_2\\
# \text{sometimes} && y(x_0) = k1 &&\text{&} && y(x_1) = k_2\\
# \end{align}
# * Let's begin by considering a simple example of a homogeneous linear $2^\circ$ ODE:
# \begin{align}
# y'' - y = 0 \implies y'' = y
# \end{align}
# Second derivative is the same as the function itself. Can you think of any solutions?
# 1. $y(x) = e^x$
# 2. $y(x) = e^{-x}$
# 3. Multiples: $y(x) = c_1 e^x$ or $y(x) = c_2e^{-x}$
# 4. Sum: $y(x) = c_1 e^x + c_2 e^{-x}$\
# check: $y' = c_1 e^x - c_2 e^{-x}$ \
# $\hspace{1cm}$ $y'' = c_1e^x + c_2 e^{-x} \implies y'' = y$ for all values of $c_1$ and $c_2$
# * Once we noticed that $y_1 = e^x$ and $y_2 = e^{-x}$ are solutions, it follows that the general linear combination is also a solution.
#   * Since $c_1$ and $c_2$ are arbitrary, this gives a double infinity family of solutions.
#   * $y_1(x)$ and $y_2(x)$ form a basis $\rightarrow$ they are linearly independent.
#   * Let's pick out a particular member of the $\infty$ solution family that also satisfies a set of initial conditions:
#   \begin{align}
#   y(0) = 2 && \text{and} &&y'(0) = -1
#   \end{align}
#   These specify the value and slope of the function ar $x=0$\
#   IC 1: 
#   \begin{align}
#   y(0) = &2 = c_1 e^0 + c_2 e^{-0}\\
#   &2 = c_1 + c_2 \implies c_1 = 2 - c_2
#   \end{align}
#   IC 2:
#   \begin{align}
#   y'(0) = -1 &= c_1 e^0 - c_2 e^{-0}\\
#   -1 &= c_1 - c_2\\
#   c_2 &= (2-c_2) + 1\\
#   2c_2 &= 3 \\
#   \implies c_2 = \frac{3}{2}\\
#   c_1 = \frac{1}{2}
#   \end{align}
#   then, $y(x) = \frac{1}{2} e^x + \frac{3}{2}e^{-x}$

# ### Solution Basis

# Consider $y(x) = c_1 e^x + c_2 e^{-x}$\
# or more generally $y(x) = c_1y_1(x) + c_2y_2(x)$\
# $\implies$ provided that $y_1(x)$ and $y_2(x)$ are linearly independent, they form a basis for the ODE solution. Similar concept to a basis of vectors\
# $\implies$ $y_1(x)$ and $y_2(x)$ are linearly independent as long as $\frac{y_1(x)}{y_2(x)}\neq $ constant

# Let's consider a general homogeneous linear $2^\circ$ ODE with constant coefficients
# \begin{align}
# y'' + a y' + by = 0
# \end{align}
# with arbitrary, real coefficients (not functions of x)
# * Based on our exprience with simple case, let us seek exponential solutions. Let's suppose $y(x) = e^{\lambda x}$ where $r$ is to be determined.\
# then, $y' = \lambda e^{\lambda x}$ and $y'' = \lambda^2 e^{\lambda x}$
# * Substituting these into our $2^\circ$ ODE gives:
# \begin{align}
# \lambda^2 e^{\lambda x} + a \lambda e^{\lambda x} + b e^{\lambda x} = 0\\
# (\lambda^2 + a\lambda + b) e^{\lambda x} = 0
# \end{align}
# Since $e^{\lambda x} \neq 0$ :
# \begin{align}
# \lambda^2 + a\lambda + b = 0 && \rightarrow \text{characteristic equation of the linear ODE}
# \end{align}
# Any $\lambda$ that is a root will give an ODE solution
# \begin{align}
# \lambda_1 = \frac{1}{2}(-a + \sqrt{a^2 - 4b})\\
# \lambda_2 = \frac{1}{2}(-a - \sqrt{a^2 - 4b})
# \end{align}
# * The characteristic equation is a quadratic equation with real coefficients. We will have three variations on the two roots:
# 1. real and different $(a^2 - 4b > 0)$
# 2. real but repeated $(a^2 - 4b = 0)$
# 3. complex conjugates $(a^2 - 4b < 0)$
# 

# ### Example (Case 1) Two real and different roots
# 
#  $y'' + 2y' - 3y = 0$ ; $y(0) = 1$ and $y'(0) = 4$
# 
#  $y(x) = c_1 e^{\lambda_1 x} + c_2 e^{\lambda_2 x}$ is the general solution\
# 

# Either solve characteristic equation: 
# \begin{align}
# \lambda^2 + a\lambda + b = 0\\
# \lambda^2 + 2\lambda - 3 = 0\\
# (\lambda + 3)(\lambda - 1) = 0 && \implies \lambda_1 = 1; \lambda_2 = -3
# \end{align}
# or with formula:
# \begin{align}
# \lambda &= \frac{1}{2}(-a \pm \sqrt{a^2 - 4b})\\
# &= \frac{1}{2}(-2 \pm \sqrt{4 + 12})\\
# &= -1 \pm 2 \implies \lambda_1 = 1 ; \lambda_2 = -3
# \end{align}
# then, 
# \begin{align}
# y(x) = c_1e^x + c_2 e^{-3x} && \text{(general solution)}
# \end{align}
# Let's double check that this works:
# \begin{align}
# y' = c_1e^x - 3c_2e^{-3x}\\
# y'' = c_1e^x + 9c_2e^{-3x}
# \end{align}
# ODE: 
# \begin{align}
# y'' + 2y' - 3y = 0\\
# \end{align}
# \begin{align}
# (c_1e^x + 9c_2e^{-3x}) + 2(c_1e^x - 3c_2e^{-3x}) - 3(c_1e^x + c_2e^{-3x}) = 0\\
# (1+2-3)c_1e^x + (9-6-3)c_3e^{-3x} = 0
# \end{align}
# Both coefficients on LHS are 0.\
# Now we need a particular solution
# \begin{align}
# y(0) = 1 &= c_1e^0 + c_2e^{-3\cdot 0}\\
# 1 &= c_1 + c_2 \implies c_2 = 1 - c_1\\
# y'(0) = 4 &= c_1e^0 - 3c_2e^{-3\cdot 0}\\
# 4 &= c_1 - 3(1-c_1)\\
# 7 &= 4c_1\\
# c_1 &= \frac{7}{4} \\
# c_2 &= -\frac{3}{4}\\
# \therefore y(x) &= \frac{7}{4}e^x - \frac{3}{4}e^{-3x} && \text{particular solution}
# \end{align}

# ### Example (Case 2): Real, double roots
# \begin{align}
# y'' + a y' + by= 0
# \end{align}
# 

# in the special circumstance that $b = \frac{1}{4}a^2$ the characteristic equation is $\lambda^2 + a\lambda + b = 0$ where
# \begin{align}
# \lambda &= \frac{1}{2}(-a \pm \sqrt{a^2 - 4b})\\
# &= \frac{1}{2}(-a \pm \sqrt{a^2 - 4\cdot \frac{1}{4}a^2})\\
# \lambda &= -\frac{1}{2}a \implies y(x) = c_1e^{-\frac{1}{2}ax}
# \end{align}
# That's great but it's only one solution (half of our basis). We need another linearly independent solution.
# Fortunately, there is a method to obtain a basis if one solution is known:

# ### Reduction of Order

# * To obtain $y_2(x)$ we set
# \begin{align}
# y_2(x) = u(x)y_1(x) \implies \text{need to determine $u(x)$}
# \end{align}
# then, 
# \begin{align}
# y_2' = u' \cdot y_1 + u \cdot y_1' && \text{(chain rule / product differentiation)}
# \end{align}
# and 
# \begin{align}
# y_2'' &= u'' \cdot y_1 + u' \cdot y_1' + u'\cdot y_1' + u \cdot y_1''\\
# &= u''y_1 + 2u'y_1' + uy_1''
# \end{align}
# * OK, let's plug those terms back into our general linear homogeneous $2^\circ$ ODE:
# \begin{align}
# y'' + p(x)y' + q(x)y = 0
# \end{align}
# \begin{align}
# u''y_1 + 2u'y_1' + uy_1'' + p(u'y_1 + uy_1') + q\cdot uy_1 = 0
# \end{align}
# * Rearrange strategically:
# \begin{align}
# u''y_1 + u'(2y_1' + py_1) + u(y_1'' + py_1' + qy_1) = 0
# \end{align}
# \begin{align}
# u'' + u'\left(\frac{2y_1' + py_1}{y_1}\right) = 0
# \end{align}

# * Now since we have $u''$ and $u'$ (but no $u$), we can reduce the order.\
# Let $U = u'$ and $U' = u''$. Now,
# \begin{align}
# U' + \left(\frac{2y_1'}{y_1} + p \right) U = 0 && \implies \text{separable!}\\
# \int \frac{dU}{U} = \int - \left(\frac{2y_1'}{y_1} + p \right)dx\\
# \ln U = -2 \ln y_1 - \int p dx\\
# U = \frac{1}{y_1^2} \exp \left[-\int p dx \right] 
# \end{align}
# Since $U = u'$, then
# \begin{align}
# u = \int U dx = \int \frac{1}{y_1^2}\exp \left[ -\int p dx \right] dx
# \end{align}
# Integral looks ugly but it's known since we know $p$ and already found $y_1$
# * Now that we have $u(x)$, we obtain our second solution, allowing us to complete our basis:
# \begin{align}
# y_2(x) = u(x) \cdot y_1(x)
# \end{align}
# and
# \begin{align}
# y(x) = c_1y_1(x) + c_2u(x)y_1(x)
# \end{align}
# $\rightarrow$ Now let's use reduction of order to solve a linear homogeneous $2^\circ$ ODE with double roots:\
# 

# ### Example
# 
# $y'' - y' + 0.25 y = 0$ ; $y(0) = 2$ and $y'(0) = \frac{1}{3}$
# 

# We can either
# 1. Solve characteristic equation
# \begin{align}
# \lambda^2 + a \lambda + b &= 0\\
# \lambda^2 - \lambda + \frac{1}{4} &= 0\\
# (\lambda - \frac{1}{2})^2 &= 0 \implies \lambda_1 = \lambda_2 = \frac{1}{2} 
# \end{align}
# 2. Use formula
# \begin{align}
# \lambda &= \frac{1}{2}(-a \pm \sqrt{a^2 - 4b})\\
# &= \frac{1}{2}(1 \pm \sqrt{1-1})\\
# \lambda &= \frac{1}{2}\\
# y_1(x) &= e^{\frac{1}{2}x}
# \end{align}
# need reduction of variables to find $y_2(x) = u(x)y_1(x)$
# 
# \begin{align}
# u(x) &= \int \frac{1}{y_1^2} \exp\left[ -\int p dx \right]dx\\
# &= \int \frac{1}{e^x}\exp\left[-\int -1 dx \right] dx\\
# &= \int e^{-x} \cdot e^x dx\\
# u(x) &= x
# \end{align}
# 
# \begin{align}
# \therefore y_2(x) = x e^{\frac{1}{2}x}
# \end{align}
# and 
# \begin{align}
# y(x) = c_1 e^{\frac{1}{2}x} + c_2 x e^{\frac{1}{2}x} && \text{general solution}
# \end{align}
# Now, 
# \begin{align}
# y(0) = 2 = c_1 e^0 + c_2 \cdot 0 \cdot e^0 \implies c_1 = 2\\
# y'(0) = \frac{1}{3} = \frac{1}{2}c_1 e^0 + c_2( e^0 + 0 \cdot e^0)\\
# \frac{1}{3} = \frac{1}{2} \cdot 2 + c_2 \implies c_2 = -\frac{2}{3} 
# \end{align}
# and 
# \begin{align}
# y(x) = 2 e^{\frac{1}{2}x} - \frac{2}{3} x e^{\frac{1}{2}x} && \text{particular solution}
# \end{align}
# 

# ### Example: Case 3 (Complex roots)
# 
# $y'' + ay' + by = 0$ (second order, homogeneous, linear with constant coefficients)

# 
# $\rightarrow$ complex roots when the characteristic equation $\lambda^2 + a \lambda + b = 0$ yields
# \begin{align}
# \lambda = \frac{1}{2}(-a \pm \sqrt{a^2 - 4b}) && \text{where } a^2 - 4b < 0
# \end{align}
# Then,
# \begin{align}
# \sqrt{a^2 - 4b} &= \sqrt{-1} \cdot \sqrt{4b - a^2}\\
# &= i \cdot \sqrt{4} \cdot \sqrt{b - \frac{1}{4}a^2}\\
# &= 2i\sqrt{b - \frac{1}{4}a^2}\\
# \therefore \lambda &= \frac{1}{2}\left(-a \pm 2i \sqrt{b - \frac{1}{4}a^2}\right)\\
# &= -\frac{1}{2}a \pm \frac{1}{2} \cdot 2i \sqrt{b - \frac{1}{4}a^2}\\
# \lambda &= -\frac{1}{2}a \pm i \omega && \text{where } \omega = \sqrt{b - \frac{1}{4}a^2} 
# \end{align}
# $\omega$ is always real.
# \begin{align}
# y(x) &= c_1 e^{\lambda_1 x} + c_2 e^{\lambda_2 x}\\
# &= c_1 e^{-\frac{a}{2}x} e^{i\omega x} + c_2 e^{-\frac{a}{2}x} e^{-i \omega x}
# \end{align}
# How to evaluate imaginary exponentials?
# 

# ### (Aside) Complex exponentials

# * For any complex number $z = s + it$,\
# $e^z = e^{s + it} = e^se^{it} \equiv e^s(\cos t + i\sin t)$\\
# \begin{align}
# e^{it} \equiv \cos (t) + i\sin (t) && \text{Euler formula}
# \end{align}
# * Look at $e^{it}$:
# \begin{align}
# t = 0 \implies e^{i0} = \cos (0) + i \sin(0) = 1\\
# t = \frac{\pi}{2} \implies e^{i \frac{\pi}{2}} = \cos (\frac{\pi}{2}) + i \sin(\frac{\pi}{2}) = i\\
# t = \pi \implies e^{i\pi} = \cos(\pi) + i \sin (\pi) = -1
# \end{align}
# $\rightarrow e^{it}$ oscillates back and forth into imaginary space.
# * Derivative:
# \begin{align}
# \frac{d}{dt}e^{it} &= \frac{d}{dt}(\cos t + i\sin t)\\
# &= -\sin t + i \cos t\\
# &= i (i \sin t + \cos t)\\
# \frac{d}{dt}e^{it} &= i e^{it}
# \end{align}
# So,
# \begin{align}
# y(x) = e^{-\frac{a}{2}x}(c_1 e^{i\omega x} + c_2 e^{-i \omega x})
# \end{align}
# and 
# \begin{align}
# e^{i\omega x} &= \cos(\omega x) + i\sin(\omega x)\\
# e^{-i\omega x} &= \cos(-\omega x) + i\sin(-\omega x)\\
# &= \cos(\omega x) - i \sin(\omega x) && \text{(properties of sin and cos)}
# \end{align}

# Now,
# \begin{align}
# y(x) &= e^{-\frac{a}{2}x}[c_1(\cos(\omega x) + i\sin(\omega x)) + c_2(\cos(\omega x) - i\sin(\omega x))]\\
# &= e^{-\frac{a}{2}x} [(c_1 + c_2) \cos(\omega x) + (c_1 - c_2) i \sin(\omega x)]
# \end{align}
# Looks messy, but remember, we can pick any basis as long as it solves the ODE.
# * If
# \begin{align}
# y_1 = e^{-\frac{a}{2}x}e^{i\omega x} = e^{-\frac{a}{2}x}[\cos(\omega x) + i\sin(\omega x)]\\
# y_2 = e^{-\frac{a}{2}x}e^{-i\omega x} = e^{-\frac{a}{2}x}[\cos(\omega x) - i\sin(\omega x)]
# \end{align}
# Let's define 
# \begin{align}
# y_3 &= \frac{1}{2}(y_1 + y_2) \implies\text{OK because it's a linear combo}\\
# &= \frac{1}{2}e^{-\frac{a}{2}x}[\cos\omega x + i\sin\omega x + \cos\omega x -i\sin\omega x]\\
# &= \frac{1}{2}e^{-\frac{a}{2}x} \cdot 2 \cos\omega x\\
# y_3 &= e^{-\frac{a}{2}x} \cdot \cos\omega x \implies \text{real solution}
# \end{align}
# and 
# \begin{align}
# y_4 &= -\frac{1}{2}i(y_1 - y_2)\\
# &= -\frac{1}{2}ie^{-\frac{a}{2}x}[\cos\omega x + i\sin\omega x - \cos\omega x + i\sin\omega x]\\
# &= -\frac{1}{2} i e^{-\frac{a}{2}x}\cdot 2i\sin\omega x\\
# y_4 &= e^{-\frac{a}{2}x} \sin\omega x \implies \text{real solution}
# \end{align}
# Then, we can say
# \begin{align}
# y(x) &= Ay_3(x) + B y_4(x) && \text{general solution}\\
# y(x) &= e^{-\frac{a}{2}x}[A\sin(\omega x) + B\cos(\omega x)] && \omega = \sqrt{b - \frac{1}{4}a^2}
# \end{align}
# where $A$ and $B$ are our arbitrary constants. All components of the solution are real.

# In[ ]:




