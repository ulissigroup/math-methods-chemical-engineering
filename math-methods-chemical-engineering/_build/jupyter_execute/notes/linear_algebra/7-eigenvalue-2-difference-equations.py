#!/usr/bin/env python
# coding: utf-8

# {{ badge }}
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$

# # $\underline{Ex}$: Symmetric matrix ($\lambda$ always real)

# * A real square matrix $\arr{A}=[a_{jk}]$ is symmetric if $\arr{A}^T = \arr{A}\implies$ thus $a_{kj}=a_{jk}$\
# $\arr{A}=\begin{bmatrix} -5&2\\2&-2 \end{bmatrix}$\
# Eigenvalues:
# \begin{align}
# |\arr{A}-\lambda\arr{I}|=0\\
# \left|\begin{array}{} -5-\lambda&2 \\2&-2-\lambda \end{array} \right|=0\\
# &=(-5-\lambda)(-2-\lambda)-4\\
# &=10+5\lambda+2\lambda+\lambda^2-4\\
# &=\lambda^2+7\lambda+6\\
# &=(\lambda+6)(\lambda+1)\\
# \implies \lambda_1=-6 ;\hspace{1cm} \lambda_2 = -1
# \end{align}
# * The eigenvalues of symmetric matrices are always real.\
# Eigenvectors:\
# For $\lambda_1=-6$, 
# \begin{align}
# (\arr{A} + 6\arr{I})\vec{x}^{(1)} = \vec{0}\\
# \begin{bmatrix} 1&2\\2&4 \end{bmatrix} \rightarrow \left[ \begin{array}{rr|r} 1&2&0\\0&0&0\end{array} \right]\\
# \implies x_1 + 2x_2 = 0\\
# x_1 = -2x_2\\
# x^{(1)}=[-2, 1]^T 
# \end{align}
# 
# For $\lambda_2=-1$, 
# \begin{align}
# (\arr{A} + \arr{I})\vec{x}^{(2)} = \vec{0}\\
# \begin{bmatrix} -4&2\\2&-1 \end{bmatrix} \rightarrow \left[ \begin{array}{rr|r} 0&0&0\\2&-1&0\end{array} \right]\\
# \implies 2x_1 - x_2 = 0\\
# x_2 = 2x_1\\
# x^{(1)}=[1, 2]^T 
# \end{align}
# Eigenspace:
# \begin{align}
# \vec{x}=\vec{0}; && \lambda=-6, \vec{x} = \begin{bmatrix} -2\\1 \end{bmatrix}; && \lambda=-1, \vec{x} = \begin{bmatrix} 1\\2 \end{bmatrix}
# \end{align}

# 
# # $\underline{Ex}:$ A skew symmetric matrix ($\lambda=0$ or complex)
# $\rightarrow \arr{A}^T = -\arr{A}$  $\hspace{0.5cm}(a_{ij}=-a_{ji})$\
# $\arr{A} = \begin{bmatrix} 0&9&-12 \\ -9&0&20 \\ 12&-20&0 \end{bmatrix}$
# * Find eigenvalues by solving $|\arr{A} - \lambda \arr{I}| = 0$
# \begin{align}
# \left|\begin{array}{} -\lambda&9&-12 \\ -9&-\lambda&20 \\ 12&-20&-\lambda \end{array}
# \right| = 0\\
# -\lambda \left|\begin{array}{} -\lambda&20\\-20&-\lambda \end{array}\right|
# -9 \left|\begin{array}{} -9&20\\12&-\lambda \end{array}\right|
# -12 \left|\begin{array}{} -9&-\lambda\\12&-20 \end{array}\right| = 0\\
# -\lambda(\lambda^2+400) -9 (9\lambda - 240) -12 (180 + 12\lambda) = 0\\
# -\lambda^3 - 400\lambda - 81\lambda + 2160 - 2160 - 144\lambda = 0\\
# -\lambda^3 - 625\lambda = 0\\
# \lambda(\lambda^2 + 625) = 0\\
# \end{align}
# \begin{align}
# \lambda_1= 0 && \lambda_2 = 25i && 
# \lambda_2 = -25i
# \end{align}
# * The eigenvalues of skew-symmetric matrices are always complex or zero.
# * Find eigenvectors:
#   1. Find $\vec{x}^{(1)}$ from $(\arr{A} - 0\arr{I}) = \vec{0}$
#   \begin{align}
#   \left[ \begin{array}{rrr|r} 0&9&-12&0 \\ -9&0&20&0 \\ 12&-20&0&0 \end{array}
#   \right]
#   \end{align}
#   Swap rows:
#   \begin{bmatrix}
#   12&-20&0 \\ -9&0&20 \\ 0&9&-12
#   \end{bmatrix}
#   $R_3 = \frac{1}{3}R_3$ ; $R_1 = \frac{1}{4}R_1$:
#   \begin{bmatrix}
#   3&-5&0 \\ -9&0&20 \\ 0&3&-4
#   \end{bmatrix}
#   $R_2 = R_2 +3R_1$, $R_3 = 5R_3$:
#   \begin{bmatrix}
#   3&-5&0 \\ 0&-15&20 \\ 0&15&-20
#   \end{bmatrix}
#   $R_3 = R_3 + R_2$:
#   \begin{bmatrix}
#   3&-5&0 \\ 0&-15&20 \\ 0&0&0
#   \end{bmatrix}
#   $R_2 = \frac{1}{5} R_2$\
#   \begin{align}
#   \left[ \begin{array}{rrr|r} 3&-5&0&0 \\ 0&-3&4&0 \\ 0&0&0&0 \end{array}
#   \right]
#   \end{align}
#   \begin{align}
#   \implies 3x_1 - 5x_2 = 0\\
#   -3x_2 + 4x_3 = 0\\
#   \implies x_2 = \frac{4}{3}x_3\\
#   x_1 = \frac{20}{9}x_3\\
#   \end{align}
# 
#   \begin{align}
#   \lambda_1 = 0, && \arr{x}^{(1)} = \begin{bmatrix}20\\12\\9 \end{bmatrix}
#   \end{align}
#   Note these are the off-diagonal terms\
#   $\vec{x}^{(2)}$ and $\vec{x}^{(3)}$ will be complex. Practice finding them @ home.

# * Eigenspace of $\arr{A}$ is:
# \begin{align}
# \vec{x} = \vec{0}; && 0, \begin{bmatrix} 20 \\ 12 \\9 \end{bmatrix}; && 25i, \begin{bmatrix} -0.33 + 0.55i \\ -0.20 - 0.92i \\ 1 \end{bmatrix} && -25i, \begin{bmatrix} -0.33 - 0.55i \\ -0.20 + 0.92i \\ 1\end{bmatrix} 
# \end{align}
# Note the complex conjugate vectors.
# 

# # Final Example with Triangular Matrix:
# 
# 

# $\arr{A} = \begin{bmatrix} 1&0&0 \\ -9&2&0 \\ 12&1&-3 \end{bmatrix} \rightarrow$ lower triangular matrix, $a_{ij} = 0$ if $j>1$
# * Eigenvalues from $|\arr{A}-\lambda \arr{I}| = 0$
# \begin{align}
# \left|
# \begin{array}{}
# 1-\lambda&0&0 \\ -9&2-\lambda&0 \\ 12&1& -3-\lambda
# \end{array}
# \right| = 0\\
# (1-\lambda)\left| \begin{array}{} 2-\lambda&0\\1&3-\lambda \end{array}\right|=0\\
# (1-\lambda)(2-\lambda)(3-\lambda)=0\\
# \end{align}
# \begin{align}
# \implies \lambda_1 = 1,&& \lambda_2 = 2,&& \lambda_3 = -3
# \end{align}
# * For upper or lower triangular matrices, eigenvalues will be diagonal elements.

# # One Application of the Eigenvalue Problem

# * "Difference" Equations (Recurrance Relations)
#   * Can be used to solve problems in population dynamics, digital signal processing and economics.
#   * We will use one here to determine how quickly people die from a plague.
# * First we must understand that the eigenvectors of any n-dimentional problem form a basis for $R^n$. In 3D space for example, we often think of $[1, 0, 0], [0, 1, 0]$ & $[0, 0, 1]$ as basis vectors of $R^n$. Eigenvectors also form a basis because they are always linearly independent.
# 
# <img src="https://docs.google.com/drawings/d/e/2PACX-1vQOiIlP4wyJlZRcYfFNs3iyuMo8zzmDVS5fe3DFiW2ZtPgw01dv-t7nYIslQ9CtOnSt8wTzmBQmBUZH/pub?w=960&h=720">
# 
# * Because eigenvectors form a basis, we can write any vector, $\vec{z}$, as 
# \begin{align}
# \vec{z} = c_1 \vec{x_1}+ c_2 \vec{x_2} + ... c_n \vec{x_n} && \text{for scalars $c_1 ... c_n$}\\
# \end{align}
# 
# \begin{align}
# e.g. [2,4,3] &= 2[1, 0, 0] + 4[0, 1, 0] + 3[0, 0, 1]\\
# &=-25[-2, -1, -1] -13[3, 3, 1] -9[1, -2, 1]
# \end{align}
# 
# \begin{align}
# \therefore \arr{A}\vec{z} &= \arr{A}(c_1 \vec{x_1} + c_2 \vec{x_2}+...c_n \vec{x_n})\\
# &= c_1\arr{A}\vec{x_1} + c_2\arr{A}\vec{x_2} + ... c_n\arr{A}\vec{x_n}
# \end{align}
# 
# 

# But, since $\arr{A}\vec{x_n} = \lambda_n\vec{x_n}$, we can greatly simplify:
# \begin{align}
# \arr{A}\vec{z} = c_1\lambda_1\vec{x_1} + c_2\lambda_2\vec{x_2} + ... + c_n\lambda_n \vec{x_n}
# \end{align}
# $\therefore$ If we know $\arr{A}$, finding its eigenvalues and vectors will help us find $\vec{z}$
# * This becomes interesting when wanting to multiply an eigenvector $\vec{x}$ by powers of $\arr{A}$. For example:
# \begin{align}
# \arr{A}^2\vec{x} = \lambda^2\vec{x} && (\arr{A}\arr{A}\vec{x} = \arr{A}\lambda\vec{x} = \lambda\arr{A}\vec{x} = \lambda\lambda\vec{x} = \lambda^2\vec{x})\\
# \arr{A}^3\vec{x} = \lambda^3\vec{x}\\
# \arr{A}^k\vec{x} = \lambda^k\vec{x}
# \end{align}
# then, $\arr{A}^k\vec{z} = c_1\lambda_1^k\arr{x_1} + c_2\lambda_2^k\vec{x_2} + ... + c_n\lambda_n^k\vec{x_n}$
# * We are interested in equations of the form $\arr{A}^k\vec{z}$ because it can tell us how some initial state represented by vector $\vec{z}$ changes on a recurrant basis. For example, if $\arr{z_0}$ represents the initial state of a population and $\arr{A}$ tells us how the population changes each year then after 1 year: 
# \begin{align}
# \vec{z_1} = \arr{A}\vec{z_0} 
# \end{align}
# And after 2 years:
# \begin{align}
# \vec{z_2} &= \arr{A}\vec{z_1} = \arr{A}\arr{A}\vec{z_0}\\
# \vec{z_2} &= \arr{A}^2\vec{z_0}
# \end{align}
# And after k years:
# \begin{align}
# \vec{z_k} = \arr{A}^k\vec{z_0} &&\rightarrow \text{this type of equation is called a difference equation} 
# \end{align}

# $\textbf{Problem:}$ A fearsome new strain of Ebola strikes the island of Niihau, HI with population of 200. People can be categorized as healthy, sick or dead. Each year, the plague causes 60% of healthy people to get sick and another 10% of healthy people to die. Only 30% stay healthy. This strain of Ebola is difficult to cure, and so 60% of sick people die each year, 20% become healthy and 20% will remain sick. 
#   * Determine the equation that predicts how many people are healthy, sick and dead after k years.
#   * How many years will it take before only 20 people remain alive on Niihau?

# * Begin by setting up variables and equations:\
# Let $z_1(k)$ = number of healthy people after k years\
# $z_2(k)$ = number of sick people after k years\
# $z_3(k)$ = number of dead people after k years\
# $\implies$ 3 types of people, 3 dimensions
# 
# * Then info in problem tells us:\
# HEALTHY
# \begin{align}
# z_1(k+1) = 0.3z_1(k) + 0.2 z_2(k)
# \end{align}
#   * In a given year $(k+1)$, the number of healthy people equals 30% of the healthy people from the previous year $(k)$ + 20% of the sick people from the previous year $(k)$.
# 
#   SICK
# \begin{align}
# z_2(k+1) = 0.6 z_1(k) + 0.2 z_2(k)
# \end{align}
# 
#   DEAD
# \begin{align}
# z_3(k+1) = 0.1 z_1(k) + 0.6 z_2(k) + z_3(k)\\ & \rightarrow \text{All dead people stay dead}
# \end{align}

# So if $\vec{z_k} = \begin{bmatrix} z_1(k) \\ z_2(k) \\ z_3(k) \end{bmatrix}$, then this problem is asking us to solve the difference equaiton 
# \begin{align}
# \vec{z_k} = \arr{A}\vec{z_{k-1}} = \arr{A}^k\vec{z_0}
# \end{align}
#   * Let's write the initial population vector in terms of eigenvalues and vectors of $\arr{A}$.
# \begin{align}
# \vec{z_k} = c_1\lambda_1^k\vec{x_1} + c_2\lambda_2^k\vec{x_2} + c_3\lambda_3^k\vec{x_3}
# \end{align}
# 
#   * Need to find eigenvalues and eigenvectors of $\arr{A} = \begin{bmatrix} 0.3&0.2&0 \\ 0.6&0.2&0 \\ 0.1&0.6&1 \end{bmatrix}$
#   \begin{align}
#   \lambda_1 = -0.1, && \vec{x}^T = [1, -2, 1]\\
#   \lambda_2 = 0.6, && \vec{x}^T = [-2, -3, 5]\\
#   \lambda_1 = 1, && \vec{x}^T = [0, 0, 1]
#   \end{align}
#   * Finally we can find the $c_n$ values by expressing the initial population vector, $\vec{z_0} = \begin{bmatrix} 200\\0\\0 \end{bmatrix}$ in terms of the eigenvector basis.
# \begin{align}
# \vec{z_0} &= \frac{600}{7}\vec{x_1} - \frac{400}{7}\vec{x_2} + 200 \vec{x_3}\\
# \begin{bmatrix}200\\0\\0\end{bmatrix} &= \frac{600}{7}\begin{bmatrix}1\\-2\\1\end{bmatrix} - \frac{400}{7}\begin{bmatrix}-2\\-3\\5\end{bmatrix} + 200 \begin{bmatrix}0\\0\\1\end{bmatrix}
# \end{align}
# Here the coefficients represent $c_1, c_2$ and $c_3$ respectively.
# \begin{align}
# \therefore \vec{z_k} = \frac{600}{7}(-0.1)^k\begin{bmatrix}1\\-2\\1\end{bmatrix} - \frac{400}{7}(0.6)^k\begin{bmatrix}-2\\-3\\5\end{bmatrix} + 200(1)^k\begin{bmatrix}0\\0\\1\end{bmatrix}
# \end{align}

#    * Finally, to determine how long it will take to have only 20 people remaining, we can plug in some values of $k$.
#    * Rounded to the nearest person:
#   \begin{align}
#    \vec{z_0} = \begin{bmatrix} 200\\0\\0 \end{bmatrix} && \vec{z_1} = \begin{bmatrix} 60\\120\\20 \end{bmatrix}\\
#    \vec{z_2} = \begin{bmatrix} 42\\60\\98 \end{bmatrix} && \vec{z_3} = \begin{bmatrix} 25\\37\\138 \end{bmatrix}\\
#    \vec{z_4} = \begin{bmatrix} 15\\22\\163 \end{bmatrix} && \vec{z_5} = \begin{bmatrix} 9\\13\\178 \end{bmatrix}\\
#    \vec{z_6} = \begin{bmatrix} 5\\8\\187 \end{bmatrix}
#   \end{align}
#   * After 6 years, only 13 people remain (less than 10% of original population)
# 
#   <img src="https://docs.google.com/drawings/d/e/2PACX-1vSjI8TG1ArprwlXWkPPKSoIacb_CbaemjDjFTcYCBfr-YdSEPhM02mitp1NW8uwvcDKh_xDrT_24DdS/pub?w=960&h=720">
