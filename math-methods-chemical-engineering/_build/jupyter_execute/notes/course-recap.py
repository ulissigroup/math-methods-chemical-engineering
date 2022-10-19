#!/usr/bin/env python
# coding: utf-8

# {{ badge }}
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$
# 

# # 06-262 Math Methods Recap
# 
# * Taught you methods needed to solve steady-state and non-steady-state problems.
# * Concepts built on each other. Matrix algebra and the concept of eigenvalues and eigenvectors became important in solving systems of ODEs. Our understanding of $1^\circ$ ODEs helped us solve $2^\circ$ ODEs. Our understanding of ODEs helps solve PDEs. And homogeneous solutions help us solve non-homogeneous problems

# ## Linear Algebra Recap
# 

# \begin{align}
# \arr{A} = \begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{bmatrix} \rightarrow 2\times 2 \text{ matrix}
# \end{align}
# * Allows us to use shorthand in defining linear systems of equations
# \begin{align}
# a_{11}x_1 + a_{12}x_2 + ... + a_{1n}x_n = b_1\\
# a_{11}x_1 + a_{12}x_2 + ... + a_{1n}x_n = b_1\\
# .\\.\\.\\
# a_{m1}x_1 + a_{m2}x_2 + ... + a_{mn}x_n = b_m\\
# \implies \arr{A}\vec{x} = \vec{b}
# \end{align}
# Which is homogneous when $\vec{b} \equiv \vec{0}$ and non-homogeneous when $\vec{b} \neq \vec{0}$
# * We used Gauss Elimination to solve systems of equations.\
# E.g.the system
# \begin{align}
# x_1 + 2x_2 - x_3 = 4\\
# 4x_2 - 2x_3 = -2\\
# x_1 - 2x_2 + 3x_3 = 0
# \end{align}
# can be represented by
# \begin{align}
# \left[\begin{array}{rrr|r} 1 & 2 & -1 & 4\\
# 0 & 4 & -2 & -2 \\
# 1 & -2 & 3 & 0\end{array}\right] \rightarrow ^{\text{G.E.}}
# \left[\begin{array}{rrr|r} 1 & 0 & 0 & 5\\
# 0 & 1 & 0 & -2 \\
# 0 & 0 & 1 & 3\end{array}\right]\\
# \implies x_1 = 5, \ x_2 = -2, \ x_3 = 3
# \end{align}

# * Systems of equations can be
#   1. underdetermined (m<n)
#   2. determined (m=n)
#   3. overdetermined (m>n)\
# And can have 3 solution scenarios:
#   1. infinite
#   2. unique
#   3. none\
# Not dictated by how "determined" system is
# 
# * Determinants $\rightarrow$ important property to calculate
# \begin{equation}
#   |\arr{A}| = \det(\arr{A})=\begin{cases}
#     a_{11} & \text{if $n=1$}.\\
#     a_{11}A_{11} + a_{12}A_{12} + ... + a_{1n}A_{1n} & \text{if $n\neq 1$}.
#   \end{cases}
# \end{equation}
# (only for square matrices) Here $A_{11}, \ A_{12}, \ A_{1n} ...$ are cofactors\
# e.g. for a $2\times 2$ matrix, $\det(\arr{A}) = a_{11}a_{22} - a_{21}a_{12}$
# 
# * A matrix is singular (non-invertible) when the det. is zero

# * Eigenvalues and eigenvectors\
# For an $n\times n$ matrix $\arr{A}$, when $
# \arr{A}\vec{x}=\lambda\vec{x}$, $\lambda$ is an eigenvalue and $\vec{x}$ is the corresponding eigenvector
#   * $\lambda$'s found by solving $|\arr{A}-\lambda\arr{I}|=0$
#   * corresponding $\vec{x}$ found by solving $(\arr{A}-\lambda\arr{I})\vec{x}=0$
#   * the eigenspace of $\arr{A}$ is $\vec{x}=\vec{0}$ (always a solution) + all $\lambda,\vec{x}$ pairs
#   * An $n\times n$ matrix yields $n$ $\lambda,\vec{x}$ pairs.
# 
# * Special matrices
#   * Diagonal
#   * Identity 
#   * Triangular (upper/lower)
#   * Symmetric and skew-symmetric

# ## ODE Recap

# * We discussed several types of ODEs. Your ability to solve requires you to classify the ODE
#   * these 5 classifications, together, will determine your solution method
# 1. Order\
# Determined by the order of the highest derivative\
#   A. first $\ y'=5$\
#   B. second $\ y'' + y' + y = 5$
# 2. Seperable-ness
#   * applies only to $1^\circ$ ODEs
#   * can be written as $g(y)\cdot y'=f(x)$\
#   A. separable $\ y^2\cdot y'= x + 2$\
#   B. non-separable $\ y' + y = x$
# 3. Exactness
#   * applies only to $1^\circ$ ODEs
#   * can be written as $N(x,y)y'=M(x,y)$ with $\frac{\partial N}{\partial x} = \frac{\partial M}{\partial y}$\
#   A. Exact $\ y\cdot y'=x^2$ (all separable ODEs are exact)\
#   B. Non-exact $y\cdot y' = x^2 + y$
# 4. Linearity
#   * can be written as $y''+p(x)y'+q(x)y=r(x)$ where $p,q$ and $r$ are continuous functions of $x$ only\
#   * Constant coefficient is case when all p(x), q(x), etc are just constants
#   A. Linear $\ y' + xy = 0$\
#   B. Non-linear $\ y\cdot y' = x^2 + y$
# 5. Homogeneity
#   * classification applies only to linear equations
#   * occurs when $r(x) \equiv 0$\
#   A. Homogeneous $y' + xy = 0$\
#   B. Non-homogeneous $y'' + y = 5x$

# #### Your Approach
# 1. If it is first order:\
#   A. separate if possible\
#   B. check linearity. If linear, make exact with integrating factor $F(x) = \exp[\int p(x)dx]$ and use chain rule to obtain $y(x) = e^{-\int p(x)dx}[\int F(x)r(x)dx + c]$\
#   C. Write in form $M(x,y)dx+N(x,y)dy=0$, check for exactness.\
#   If not exact, try to find integrating factor $F(x)$ or $F(y)$
#   * If $F(x)$, then $\frac{1}{F}\frac{dF}{dx} = \frac{1}{N}\left\{\frac{\partial M}{\partial y} - \frac{\partial N}{\partial x}\right\} = R(x)$\
#   and $F(x)=\exp[\int R(x)dx]$
#   * If $F(y)$, then $\frac{1}{F}\frac{dF}{dy} = \frac{1}{M}\left\{\frac{\partial N}{\partial x} - \frac{\partial M}{\partial y}\right\} = R(y)$\
#   and $F(y)=\exp[\int R(y)dy]$\
#   Once your ODE is exact, find $u(x,y)$ s.t. $\frac{\partial u}{\partial x} = M(x,y)$ and $\frac{\partial u}{\partial y} = N(x,y)$\
#   Implicit solution is given by $u(x,y)=c$
# 2. If it is second order:\
#   A. check homogeneity. If homogeneous, identify the two component solution basis, which depends on the roots of the characteristic equation $\lambda^2+a\lambda+b=0$
#   1. real, distinct roots\
#   $y(x)=c_1e^{\lambda_1x}+c_2e^{\lambda_2x}$, where $\lambda=\frac{1}{2}(a\pm\sqrt{a^2-4b})$
#   2. repeated roots\
#   $y(x)=c_1e^{\lambda x}+c_2xe^{\lambda x}$ where $\lambda=-\frac{1}{2}a$
#   3. complex roots\
#   $y(x)=e^{-\frac{a}{2}x}[c_1\sin(\omega x)+c_2\cos(\omega x)]$ where $\omega=\sqrt{b-\frac{1}{4}a^2}$
#   
#   B. if non-homogeneous, find homogeneous solution $y_H(x)$, then find a particular solution $y_P(x)$ using either
#   1. Method of undetermined coefficients
#     * applicable only for non-homogeneous terms that return themselves as derivatives, e.g. $\sin x, e^x, x^n$
#     * assume $y_P(x)$ looks like $r(x)$
#   2. Variation of Parameters
#     * works for any $r(x)$
#     * assume $y_P(x)=u(x)y_1(x)+v(x)y_2(x)$
#     * Find Wronskian, $W=y_1y_2'-y_2y_1'$, then,\
#     $y(x)=y_H(x)+y_1(x)\int\frac{-r(x)y_2(x)}{W}dx + y_2(x)\int\frac{r(x)y_1(x)}{W}dx$
# 
# 
# * Boundary Value Problems
#   * IVPs always have a unique solution provided continuity is met
#   * BVPs have either no solution, unique or infinite number of solutions $\rightarrow$ you must think for these

# ## Coupled ODEs
# 
# 
# 

# e.g.\begin{align}
# y_1'-y_1+y_2=4\\
# y_2'+y_1=8x
# \end{align}
# especially when working with higher number of equations (like in chemical plant), best to solve simultaneously
# * Matrices come back:
# \begin{align}
# \vec{y}'=\arr{A}\vec{y}+\vec{b}(x)
# \end{align}
# where $\vec{b}(x)$ is the non-homogeneous term, $\vec{y}$ is an $n\times 1$ vector of unknown functions and $\arr{A}$ is an $n\times n$ coefficient matrix\
# For the example, $\arr{A}=\begin{bmatrix}1&-1 \\ -1&0 \end{bmatrix}$ and $\vec{b}=\begin{bmatrix} 4 \\ 8x \end{bmatrix}$
# * We proved that the eigenvalues and eigenvectors of $\arr{A}$ yield the homogenoeus solution:
# \begin{align}
# \vec{y} = c_1\vec{x}^{(1)}e^{\lambda_1t} + c_2\vec{x}^{(2)}e^{\lambda_2t} + ... + c_n\vec{x}^{(n)}e^{\lambda_nt}
# \end{align}
# for real, distinct eigenvalues
# * For real, repeated eigenvalues, we must use reduction of order to obtain linearly independent solutions\
# e.g. $c_1\vec{x}^{(1)}e^{\lambda_1t} + c_2\vec{x}^{(1)}te^{\lambda_1t}$ (remember there is also a term from the solution of the generalized eigenvalue problem!)
# * For complex eigenvalues, we must use Euler's formula
# \begin{align}
# e^{\omega it} = \cos(\omega t) + i\sin(\omega t)
# \end{align}
# to change basis and obtain a real solution
# \begin{align}
# \vec{y}=c_1 \begin{bmatrix} \cos\omega t \\ -sin\omega t \end{bmatrix} + c_2 \begin{bmatrix} \sin \omega t \\ \cos \omega t \end{bmatrix}
# \end{align}
# * Non-homogeneous systems are handled similarly to non-homogeneous single equations
#   * after finding homogeneous solution, find a particular solution using MUC or VoP
#   * In both cases, you assume a form of $y_P$\
#   MUC: $\vec{y}_P(t) = \vec{u}\sin t + \vec{v}\cos t$ (for, say, $5\cos t$)\
#       $\vec{y}_P(t) = \vec{u}t^2 + \vec{v}t + \vec{w}$ (for, say, $2t^2$)\
#   VoP: $\vec{y}_P(t) = \arr{Y}(t)\vec{u}(t)$\
#   where $\arr{Y}(t)$ is the fundamental matrix obtained from your two linearly independent solutions to homogeneous equation
# * For non-linear ODE's where we can't get the solution, we can still find steady states and do a non-linear stability analysis
#   * Steady states are straightforward $\vec{y}'=\vec{0}$, but you have to find **every steady state!**. Be careful. 
#   * At each steady state, evaluate the jacobian and form a linearized ODE that holds near the steady state. 
#   * Calculate the eigenvalues to understand the type (stability) of the steady state

# # PDEs
# 
#   * 2 or more independent variables
#   * Order, linearity, homogeneity
#   * solve using separation of variables to turn into multiple ODE's, or if given special solution forms for specific PDE's
#   * If needed, use Fourier coefficients to transform initial conditions into the basis of solutions
# 

# ### Statistics 
# 
# At the very end, we covered some simple statistics that will help you next year. 
# * Basic statistics
#   * Types of distributions (continuous, discrete) and how to evaluate them (PDF, CDF, sampling) in python
#   * Plotting distribution (**NOT JUST HISTOGRAMS**) and kernel density estimation. 
#   * Fitting distributions based on the CDF
#   * Estimates for properties of a Gaussian distribution, and the errors on those estimates
# * Regression
#   * Linear regression **[testable]**
#     * Form the augmented matrix $\arr{X}$ with your non-linear functions
#     * Solve for the best parameters using the $\arr{X}^T$ trick
#     * If you need uncertainty, use a package like statsmodels
#   * Non-linear regression
#     * Turn a non-linear problem into a linear one, then solve like that
#     * Use a tool like curvefit (easiest) or lmfit (if need uncertainty)
# 
