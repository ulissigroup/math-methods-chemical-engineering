#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/ulissigroup/math-methods-chemical-engineering/blob/master/lecture_notes/06-eigenvalues-1.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$

# Note: finished about on time, even with bringing Laika to class. 

# # Homogeneous vs. Non-homogeneous Systems

# System of linear equations represented as:
# \begin{align}
# \arr{A}\vec{x} = \vec{b}
# \end{align}
# * Non-homogeneous system if $\vec{b}\neq\vec{0}$
#   * If $det(\arr{A}\neq0)$ (e.g. $\arr{A}$ is invertible), then a unique solution exists.
# 
# * Homogeneous system if $\vec{b}=\vec{0}$
#   * If $det(\arr{A}\neq0) \rightarrow$ has only the trivial solution $\vec{x} = \vec{0}$
#   * If $det(\arr{A}=0) \rightarrow$ also has a series of nontrivial solutions

# # Eigenvalues and Eigenvectors

# * For engineering applications, eigenvalue problems are among the most important problems concerning matrices.
#   * For example, wherever there are vibrations, there are eigenvalues, which are the natural frequencies of the vibration.
#   * If you've ever tuned a guitar, you've solved an eigenvalue problem!
#   * Google's page-rank algorithm for determining which pages are important is also built on eigenvectors of a very specific matrix.
# 
# * For the following expression:
# \begin{align}
# \arr{A}\vec{x}=\lambda \vec{x} 
# \end{align}
# $\vec{x}=0$ is always a solution (trivial)\
# If $\vec{x}\neq0$ exists, then $\lambda=$ eigenvalue of $\arr{A}$ and $\vec{x}=$ eigenvctor of $\arr{A}$.
# * Often, there are are many eigenvectors, which together with the $\vec{0}$ vector, form the eigenspace of $\arr{A}$.
# * How do we find $\lambda$'s and associated $\vec{x}$'s? What do they mean?

# * To find $\lambda$, rearrange our key equation:
# \begin{align}
# \arr{A}\vec{x} = \lambda\vec{x}\\
# \arr{A}\vec{x} - \lambda\vec{x} = \vec{0}
# \end{align}
# $\arr{A}$ is matrix and $\lambda$ a scalar, but we can factor out $\vec{x}$ using $\arr{I}$.
# \begin{align}
# (\arr{A} - \lambda\arr{I})\vec{x}= \vec{0}
# \end{align}
# This is a homogeneous system. We just learned that if the determinant of the matrix on the LHS is zero, then there is a non-trivial solution that is actually a series.\
# $\therefore$ if $|\arr{A} - \lambda\arr{I}|=0$,\
# then there is a non trivial value of $\vec{x}$ and a set of $\lambda$, $\vec{x}$ that satisfy $\arr{A}\vec{x}=\lambda\vec{x}$
# 

# Ex: $\arr{A} = \begin{bmatrix}-5 & 0\\ 1 & 2\end{bmatrix}$\
# We want 
# \begin{align}
# |\arr{A} - \lambda\arr{I}|=0\\
# \left|
# \begin{bmatrix}
# -5 & 0\\
# 1 & 2
# \end{bmatrix}
# -\lambda 
# \begin{bmatrix}
# 1&0\\
# 0&1
# \end{bmatrix}
# \right| = 0\\
# \\
# \left|
# \begin{array}{}
# -5-\lambda&0\\
# 1&2-\lambda
# \end{array}
# \right| = 0
# \end{align}

# \begin{align}
# (-5-\lambda)(2-\lambda)-0 = 0\\
# (-5-\lambda)(2-\lambda)= \lambda^2 + 3 \lambda - 10 = 0 \\
# \end{align}
# which is the "characteristic equation" of $\arr{A}$.
# * Equation is satisfied by $\lambda_1=2$; $\lambda_2=-5$
# * Now, there is one eigenvector assiciated with each eigenvalue: $\vec{x}^{(1)}$ and $\vec{x}^{(2)}$
# * Let's start with $\lambda_1=2$:
# \begin{align}
# (\arr{A} - 2 \arr{I})\vec{x}^{(1)} &= \vec{0}\\
# \left(
# \begin{bmatrix}
# -5 & 0\\
# 1 & 2
# \end{bmatrix}
# -\begin{bmatrix}
# 2 & 0\\
# 0 & 2
# \end{bmatrix}
# \right)
# \begin{bmatrix}
# x_1^{(1)}\\
# x_2^{(1)}
# \end{bmatrix}&=
# \begin{bmatrix}
# 0\\0
# \end{bmatrix}\\
# \begin{bmatrix}
# -7 & 0\\
# 1 & 0
# \end{bmatrix}
# \begin{bmatrix}
# x_1^{(1)}\\
# x_2^{(1)}
# \end{bmatrix}&=
# \begin{bmatrix}
# 0\\0
# \end{bmatrix}\\
# -7x_1^{(1)}&=0\\
# x_1^{(1)}&=0
# \end{align}
# $\therefore x_1^{(1)} =0$ and $x_1^{(2)}$ can be anything but zero since solution is non-trivial.\
# $\therefore \vec{x}^{(1)} = \begin{bmatrix}0\\\alpha\end{bmatrix}$ where $\alpha\neq 0$. This $\vec{x}^{(1)}$ is married to $\lambda_1=2$.
# 
# 

# * Now, $\lambda_2=-5:$ 
# \begin{align}
# (\arr{A} - 5 \arr{I})\vec{x}^{(2)} &= \vec{0}\\
# \left(\begin{bmatrix}-5 & 0\\1 & 2\end{bmatrix}
# +\begin{bmatrix}5 & 0\\0 & 5\end{bmatrix}\right)
# \begin{bmatrix}x_1^{(2)}\\x_2^{(2)}\end{bmatrix}&=
# \begin{bmatrix}0\\0\end{bmatrix}\\
# \begin{bmatrix}0 & 0\\1 & 7\end{bmatrix}
# \begin{bmatrix}x_1^{(2)}\\x_2^{(2)}\end{bmatrix}&=
# \begin{bmatrix}0\\0\end{bmatrix}\\
# \end{align}
# Which gives $0=0$ and 
# \begin{align}
# x_1^{(2)} + 7 x_2^{(2)} = 0\\ 
# x_1^{(2)} = -7 x_2^{(2)}
# \end{align}
# $\therefore \vec{x}^{(2)} = \begin{bmatrix}-7\beta\\ \beta \end{bmatrix}$, where $\beta \neq 0$. This $\vec{x}^{(2)}$ is married to $\lambda_2 = -5$

# * Make sure to check tat $\lambda + \vec{x}$ 's satisfy the original equation\
# \begin{align}
# \arr{A}\vec{x} = \lambda \vec{x}\\
# \end{align}
# \begin{align}
# \lambda_1=2 && \text{and} && \lambda_2 = -5\\
# \end{align}
# \
# \begin{align}
# \begin{bmatrix}-5 & 0\\1 & 2\end{bmatrix}
# \begin{bmatrix}0\\ \alpha \end{bmatrix}=
# 2 \begin{bmatrix}0\\ \alpha \end{bmatrix} && \text{and} &&
# \begin{bmatrix}-5 & 0\\1 & 2\end{bmatrix}
# \begin{bmatrix}-7 \beta\\ \beta \end{bmatrix}=
# -5 \begin{bmatrix}-7 \beta \\ \beta \end{bmatrix}
# \end{align}
# \
# \begin{align}
# \begin{bmatrix} 0 \\ 2 \alpha \end{bmatrix} = \begin{bmatrix} 0 \\ 2 \alpha \end{bmatrix} &&\text{and}&&
# \begin{bmatrix} 35 \beta \\ -5 \end{bmatrix} = \begin{bmatrix}  35 \beta \\ -5\end{bmatrix}
# \end{align}
# 

# * Eigenspaces are:
# \begin{array}{}
# \begin{bmatrix} 0\\0 \end{bmatrix} ; & 2,\begin{bmatrix} 0\\ \alpha \end{bmatrix}; & -5, \begin{bmatrix} -7 \beta \\ \beta \end{bmatrix}
# \end{array}
# 
# The basis of two non zero eignspaces are $\begin{bmatrix} 0 \\ 1\end{bmatrix}$ and $\begin{bmatrix} -7 \\ 1\end{bmatrix}$ respectively.
# * For any n$\times$n matrix with $rank(\arr{A})=n$, you'll get an $n^{(th)}$ degree polynomial to solve in $\lambda$ from $|\arr{A} - \lambda \arr{I}| = 0$. 
# * An n$\times$n matrix has at least one $\lambda$ and $n$ distinct $\lambda$ 's.
# * Recap to solving eigenvalue problems:
#   1. Set up $|\arr{A} - \lambda \arr{I}| = 0$
#   2. Determine the characteristic equation to solve for $\lambda$
#   3. For each $\lambda$, determine $\vec{x}$
# 
# * Let's try a larger example:\
# $\arr{A} = \begin{bmatrix}6&10&6 \\ 0&8&12 \\ 0&0&2 \end{bmatrix} \rightarrow$ we want $\lambda$, $\vec{x}$ such that $\arr{A} \vec{x} = \lambda \vec{x}$
# \begin{align}
# |\arr{A} - \lambda \arr{I}|  &= \begin{bmatrix}6 - \lambda&10&6 \\ 0&8 - \lambda&12 \\ 0&0&2 - \lambda \end{bmatrix}\\
# \end{align}
# We choose to work with column 1
# \begin{align}
# &=(6-\lambda)\left|\begin{array}{} 8 - \lambda & 12 \\ 0 & 2 - \lambda \end{array} \right|
# - 0 \left|\begin{array}{} 10 & 6 \\ 0 & 2 - \lambda \end{array} \right|
# + 0 \left|\begin{array}{} 10 & 6 \\ 8 - \lambda & 12 \end{array} \right|\\
# &= (6-\lambda) [(8 - \lambda)(2-\lambda) - 0]\\
# 0 &= (6-\lambda)(8-\lambda)(2-\lambda)\\
# \end{align}
# \
# We don't really need to go further unless we want expanded characteristic equation.
# \begin{align}
# 0 &= \lambda^3 -16 \lambda^2 + 76\lambda - 96
# \end{align}
# 

# $\therefore$ eigenvalues are:
# \begin{array}{}
# \lambda_1=6, & \lambda_2=8, & \lambda_3=3
# \end{array}
# We can get maximum of 3 eigenvalues since $n=3$
# * Now, find $\vec{x}^{(i)}$ for each $\lambda_i$ using $(\arr{A}- \lambda_i \arr{I})\vec{x}^{(i)} = \vec{0}$
#   1. For $\lambda_1=6$:
#   \begin{align}
#   \begin{bmatrix} 0&10&6 \\ 0&2&12 \\ 0&0&-4 \end{bmatrix}
#   \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = 
#   \begin{bmatrix} 0 \\ 0 \\ 0 \end{bmatrix}
#   \end{align}
#   Which gives:
#   \begin{align}
#   10x_2 + 6x_3 = 0;\\
#   2 x_2 + 12 x_3 = 0;\\
#   -4x_3 = 0
#   \end{align}
#   
#   \begin{align}
#   \implies x_3 = x_2 = 0 &&\text{and}&& \text{$x_1=$arbitrary}
#   \end{align}
# 
#   \begin{align}
#   \vec{x}^{(1)} = \begin{bmatrix} \alpha \\ 0 \\ 0 \end{bmatrix} & \text{, where $\alpha\neq0$}
#   \end{align}
#   or the basis vector $\vec{x}^{(1)} = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}$ (really a series of solutions)
# 
#   2. For $\lambda_2=8$:
#   \begin{align}
#   \begin{bmatrix} -2&10&6 \\ 0&0&12 \\ 0&0&-6 \end{bmatrix}
#   \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = 
#   \begin{bmatrix} 0 \\ 0 \\ 0 \end{bmatrix}
#   \end{align}
#   Which gives:
#   \begin{align}
#   -2x_1 + 10x_2 + 6x_3= 0;\\
#   12x_3 = 0;\\
#   -6x_3 = 0
#   \end{align}
#   
#   \begin{align}
#   \implies x_3 = 0, && x_1 = 5x_2 + 3x_3 &&\text{and}&& \text{$x_2=$arbitrary}
#   \end{align}
# 
#   \begin{align}
#   \vec{x}^{(2)} = \begin{bmatrix} 5\beta \\ \beta \\ 0 \end{bmatrix}& \text{, $\beta\neq0$}
#   \end{align}
#   or the basis vector $\vec{x}^{(2)} = \begin{bmatrix} 5 \\ 1 \\ 0 \end{bmatrix}$
# 
#   3. For $\lambda_3=2$:
#   \begin{align}
#   \begin{bmatrix} 4&10&6 \\ 0&6&12 \\ 0&0&0 \end{bmatrix}
#   \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} = 
#   \begin{bmatrix} 0 \\ 0 \\ 0 \end{bmatrix}
#   \end{align}
#   Which gives:
#   \begin{align}
#   4x_1 + 10x_2 + 6x_3 = 0;\\
#   6x_2 + 12x_3 = 0;\\
#   0 = 0
#   \end{align}
#   
#   \begin{align}
#   \implies x_2 &= -2 x_3;\\
#   x_1 &= -\frac{10}{4}(-2x_3) - \frac{6}{4}x_3 \\
#   &= 5x_3 -\frac{3}{2}x_3\\
#   x_1 &= \frac{7}{2}x_3
#   \end{align}
# 
#   \begin{align}
#   \vec{x}^{(3)} = \begin{bmatrix} \frac{7}{2}\delta \\ -2\delta \\ \delta \end{bmatrix}& \text{, $\delta\neq0$}
#   \end{align}
#   or the basis vector $\vec{x}^{(3)} = \begin{bmatrix} \frac{7}{2} \\ -2 \\ 1 \end{bmatrix} = \begin{bmatrix} 7 \\ -4 \\ 2 \end{bmatrix}$

# Eigenspace corresponding to $\arr{A}$ :
# \begin{array}{}
# 6, \begin{bmatrix} 1\\0\\0 \end{bmatrix} ; & 8,\begin{bmatrix} 5\\1\\0 \end{bmatrix}; & 2, \begin{bmatrix} 7\\-4\\2\end{bmatrix} ; & \vec{x}=\vec{0}
# \end{array}

# # Example with a double root 
# 
# Consider $\arr{A}=\begin{bmatrix}-2&2&-3\\2&1&-6\\-1&-2&0\end{bmatrix}$
# 
# ## Find the characteristic equation
# 
# \begin{align*}
# |\arr{A}-\lambda\arr{I}|&=0\\
# \begin{vmatrix}
# -2-\lambda&2&-3\\2&1-\lambda&-6\\-1&-2&0-\lambda
# \end{vmatrix}&=0
# \end{align*}
# 
# Calculating this, we get $\lambda^3+\lambda^2-21\lambda-45=0$.  No way for us to calculate this easily by hand! 
# 
# Using python, we find the roots are $\lambda_1=5,\lambda_2=\lambda_3=-3$. Note the double root!
# 
# ##  Find the eigenvector for $\lambda_1=5$, corresponding to $(\arr{A}-5\arr{I})\vec{x}^{(1)}=\vec{0}$
# \begin{align*}
# \begin{bmatrix}
# -7&2&-3\\
# 2&-4&-6\\
# 1&-2&-5
# \end{bmatrix}\begin{bmatrix}x_1\\x_2\\x_3\end{bmatrix}=\begin{bmatrix}0\\0\\0\end{bmatrix}
# \end{align*}
# 
# We need to use Gauss Elimination to solve this:
# * $R_1=R_2/2, R_2=R_1$
# \begin{align*}
# \begin{bmatrix}
# 1&-2&-3\\
# -7&2&-3\\
# -1&-2&-5
# \end{bmatrix}
# \end{align*}
# * $R_2=R_2+7R_1$, $R_3=R_3+R_1$
# \begin{align*}
# \begin{bmatrix}
# 1&-2&-3\\
# 0&-12&-24\\
# 0&-4&-8
# \end{bmatrix}
# \end{align*}
# * $R_3=-R_3/4, R_2=-R_2/12$
# \begin{align*}
# \begin{bmatrix}
# 1&-2&-3\\
# 0&1&2\\
# 0&1&2
# \end{bmatrix}
# \end{align*}
# * $R_3=R_3-R_2, R_1=R_1+2R_2$
# \begin{align*}
# \begin{bmatrix}
# 1&0&1\\
# 0&1&2\\
# 0&0&0
# \end{bmatrix}
# \end{align*}
# 
# So $x_1+x_3=0$ and $x_2+2x_3=0$. Let's let $x_3=1$. That gives us $x_1=-\alpha$ and $x_2=-2\alpha$, so
# \begin{align*}
# \vec{x}^{(1)}=\begin{bmatrix}
# -1\\-2\\1
# \end{bmatrix}
# \end{align*}
# 
# We can do a quick check that this is correct: 
# \begin{align*}
# \arr{A}\vec{x}^{(1)}=\begin{bmatrix}-2&2&-3\\2&1&-6\\-1&-2&0\end{bmatrix}\begin{bmatrix}
# -1\\-2\\1
# \end{bmatrix}=\begin{bmatrix}
# -5\\-10\\5
# \end{bmatrix}=5\begin{bmatrix}
# -1\\-2\\1
# \end{bmatrix}
# \end{align*}
# 
# ## Eigenvector for $\lambda_2=\lambda_3=-3$
# 
# For these, we get:
# \begin{align*}
# \begin{bmatrix}
# 1&2&-3\\
# 2&4&-6\\
# -2&2&3
# \end{bmatrix}\begin{bmatrix}x_1\\x_2\\x_3\end{bmatrix}=\begin{bmatrix}0\\0\\0\end{bmatrix}
# \end{align*}
# After Gauss elimination we get:
# \begin{align*}
# \begin{bmatrix}
# 1&2&-3\\
# 0&0&0\\
# 0&0&0
# \end{bmatrix}\begin{bmatrix}x_1\\x_2\\x_3\end{bmatrix}=\begin{bmatrix}0\\0\\0\end{bmatrix}
# \end{align*}
# This tells us that $x_1+2x_2-3x_3=0$. This has one equation and three unknowns! So, let's write our system as a system of three equations:
# \begin{align*}
# x_1&=-2x_2+3x_2\\
# x_2&=x_2\\
# x_3&=x_3
# \end{align*}
# We can write this as three linearly independent vectors:
# \begin{align*}
# \vec{x}=\begin{bmatrix}
# x_1\\x_2\\x_3
# \end{bmatrix}=x_2\begin{bmatrix}
# -2\\1\\0
# \end{bmatrix}+x_3\begin{bmatrix}
# 3\\0\\1
# \end{bmatrix}
# \end{align*}
# where $x_2$ and $x_3$ are arbitrary. This gives us two eigenvectors:
# \begin{align*}
# \vec{x}^{(2)}=\begin{bmatrix}
# -2\\1\\0
# \end{bmatrix}&&\vec{x}^{(3)}=\begin{bmatrix}
# 3\\0\\1
# \end{bmatrix}
# \end{align*}
# 
# So, the fulll eigenspace of $\arr{A}$ is 
# \begin{align*}
# \vec{x}=\vec{0} && 5,\begin{bmatrix}
# -1\\-2\\1
# \end{bmatrix} && -3, \begin{bmatrix}
# -2\\1\\0
# \end{bmatrix},\begin{bmatrix}
# 3\\0\\1
# \end{bmatrix}
# \end{align*}

# # Recap of the characteristic equation
# 
# $|\arr{A}-\lambda \arr{I}|=0$ gives us the characteristic  equation of $\arr{A}$, which is a polynomial. An nth order polynomial has several solution possibilities:
# * $n$ distinct, real roots
# * redundant, real roots
# * complex (i.e. imaginary) roots
