#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/ulissigroup/math-methods-chemical-engineering/blob/master/lecture_notes/05-determinants.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$

# # Using the inverse of A to solve a system of linear equations

# * If you have the inverse of a matrix $\arr{A}$ representative of a linear system of equations, you can use it to find $\vec{x}$:
# \begin{align*}
# \arr{A}\vec{x}=\vec{b}\\
# \arr{A}^{-1}\arr{A}\vec{x}=\arr{A}^{-1}\vec{b}\\
# \vec{x}=\arr{A}^{-1}\vec{b}
# \end{align*}
# now $\vec{x}$ is just the product.
# * Fun fact: An n$\times$n system $\arr{A}\vec{x}=\vec{b}$ has a unique solution if and only if $\arr{A}$ is nonsingular (i.e. has an inverse).
# *Other properties:
# \begin{align*}
# (\arr{A}\arr{C})^{-1}=\arr{C}^{-1}\arr{A}^{-1}\\
# (\arr{A}\arr{B}\arr{C})^{-1}=\arr{C}^{-1}\arr{B}^{-1}\arr{A}^{-1}
# \end{align*}

# * Be careful:\
# if $\arr{A}=\begin{bmatrix}1&0&-2\\
# 2&-1&3\\
# 4&0&2\end{bmatrix}$ then $\arr{A}^{-1}=\begin{bmatrix}2/10&0&2/10\\
# -8/10&-1&7/10\\
# \bf{-4/10}&\bf{0}&\bf{1/10}\end{bmatrix}$

# if $\arr{B}=\begin{bmatrix}1&0&-2\\
# 4&0&2\\
# 2&-1&3\end{bmatrix}$(switched rows), then $\arr{B}^{-1}\neq\begin{bmatrix}2/10&0&2/10\\
# \bf{-4/10}&\bf{0}&\bf{1/10}\\
# -8/10&-1&7/10\end{bmatrix}$ but $\arr{B}^{-1}=\begin{bmatrix}2/10&2/10&0\\
# -8/10&7/10&-1\\
# -4/10&1/10&0\end{bmatrix}$
# 
# Notice the columns were switched here.
# 

# # Determinants

# * Every square matrix is associated with a real number (a scalar) called the determinant.
# * The calue of the determinant will tell us whether or not the matrix is singular.
# * Generalized definition:
# \begin{equation}
#   D = |\arr{A}| = det(\arr{A})=\begin{cases}
#     a_{11} & \text{if $n=1$}.\\
#     a_{11}A_{11} + a_{12}A_{12} + ... + a_{1n}A_{1n} & \text{if n>1}.
#   \end{cases}
# \end{equation}

# where $A_{1k} = (-1)^{1+k} det(\arr{M_{1k}})$, for $k=1,2,...n$\
# Here, $A_{1k}$ is a "cofactor" and $det(\arr{M_{1k}})$ is the "minor" of $a_{1k}$.\
# Also, $\arr{M_{1k}}$ = submatrix formed from $\arr{A}$ when deleting row $1$ and column $k$.
# 
# OK, so this makes no sense without some examples. Let's consider matrices of different sizes.

# ## Case 1: 1$\times$1 Matrix

# $|A| = a_{11}=a$
# * $\arr{A}$ is non-singular if $a\neq0$ and singular if $a=0$.\
# Ex: $\arr{A}=[5] \hspace{0.5cm} \rightarrow{} det(\arr{A})=5$ (non singular) \
# $\hspace{0.6cm} \arr{B}=[0] \hspace{0.5cm} \rightarrow{} det(\arr{B})=0$ (singular) 
# 
# Let's think about a physical system $ax=5$. If $a$ is zero, we have a problem! If it's non-zero, we can solve.

# ## Case 2: 2$\times$2 matrix

# Let $\arr{A} = \begin{bmatrix}a_{11} & a_{12}\\ a_{21} & a_{22}\end{bmatrix}$
# * $\arr{A}$ will be non-singular only if we can convert it to $\arr{I}$ using elementary row operations. This property is called "row equivalence".
# * Let's consider how we would convert $\arr{A}$ to $\arr{I}$.\
# If $a_{11}\neq0$ then we can perform the following row operations:
# 1.   Multiply the second row of $\arr{A}$ by $a_{11}$ $(R_2=a_{11}R_2)$
# 
# \begin{align*}
# \arr{A}=\begin{bmatrix}a_{11} & a_{12}\\
# a_{11}a_{21} & a_{11}a_{22}\end{bmatrix}
# \end{align*}
# 
# 2. Subtract $a_{21} \times $Row 1 from Row 2 $(R_2 = R_2 - a_{21} R_1)$\
# \begin{align*}
# \arr{A}=\begin{bmatrix}a_{11} & a_{12}\\
# 0 & a_{11}a_{22} - a_{12}a_{21}\end{bmatrix}
# \end{align*}
# 
# * Since $a_{11}\neq0$, $\arr{A}$ will be row equivalent to $\arr{I}$ IFF:\
# \begin{align*}
# a_{11}a_{22}-a_{12}a_{21} \neq 0 \hspace{2cm} (*)
# \end{align*}
# * If $a_{11}=0$, then we can switch the two rows:
# \begin{align*}
# \arr{A} = \begin{bmatrix}
# a_{21} & a_{22}\\
# 0 & a_{12}
# \end{bmatrix}
# \end{align*}
# * Now, $\arr{A}$ will be row equivalent to $\arr{I}$ if $a_{21}\neq0$ and $a_{12}\neq0$
# * This is the same as saying $a_{11}a_{21} \neq 0$ which is equivalent to (*) for $a_{11}=0$:
# \begin{align*}
# a_{11}a_{22} - a_{12}a_{21} = a_{12}a_{21}\neq0
# \end{align*}
# * Great! These cases of $a_{11}\neq0$ and $a_{11}=0$ are consistent and we can define \
# $det(\arr{A}) = |\arr{A}| = a_{11}a_{22} - a_{12} a_{21}$
# \begin{equation}
#   \begin{cases}
#     =0 & \text{when $\arr{A}$ is singular}\\
#     \neq 0 & \text{when $\arr{A}$ is non-singular}
#   \end{cases}
# \end{equation}
# 

# ## Case 3: 3$\times$ 3 Matrix

# * We will not derive the determinant for a 3 $\times$ 3 matrix in class (see text).
# * How to find it?
# \begin{align*}
# det(\arr{A}) = \left|\begin{array}{}
# a_{11}&a_{12}&a_{13}\\
# a_{21}&a_{22}&a_{23}\\
# a_{31}&a_{32}&a_{33}
# \end{array}\right|
# \end{align*}
# * Working with Row 1, decompose as follows:
# \begin{align}
# det(\arr{A}) &= a_{11} \left|\begin{array}{}
# a_{22}&a_{23}\\
# a_{32}&a_{33}
# \end{array}\right|
# -a_{12} \left|\begin{array}{}
# a_{21}&a_{23}\\
# a_{31}&a_{33}
# \end{array}\right|
# + a_{13} \left|\begin{array}{}
# a_{21}&a_{22}\\
# a_{31}&a_{33}
# \end{array}\right|\\
# &= a_{11}(a_{22}a_{33}-a_{23}a_{32}) - a_{12}(a_{21}a_{33}-a_{23}a_{31}) + a_{13}(a_{21}a_{32}-a_{22}a_{31})
# \end{align}
# * Looks complicated but it's just a number. Each determinant is a "minor" determinant of the submatrix formed by eliminating row $j$ and column $k$.\
# $Ex: \arr{A} = \left|\begin{array}{}
# 2&5&4\\
# 3&1&2\\
# 5&4&6\end{array}\right|$ then 
# \begin{align}
# det(\arr{A}) &= 2 \left|\begin{array}{}
# 1&2\\
# 4&6\end{array}\right|
#  - 5 \left|\begin{array}{}
# 3&2\\
# 5&6\end{array}\right|
# +4 \left|\begin{array}{}
# 3&1\\
# 5&4\end{array}\right|\\
# &=2(6-8) - 5(18-10) + 4(12-5)\\
# &=-4 -40 + 28\\
# det(\arr{A})&=-16
# \end{align}

# * You don't need to work with Row 1. You can use any column or row you'd like but must switch signs for even number rows or columns.\
# Ex: Choose column 2\
# $\arr{A} = \left|\begin{array}{}
# 2&5&4\\
# 3&1&2\\
# 5&4&6\end{array}\right|$ 
# \begin{align}
# det(\arr{A}) &= -5 \left|\begin{array}{}
# 3&2\\
# 5&6\end{array}\right|
#  +1 \left|\begin{array}{}
# 2&4\\
# 5&6\end{array}\right|
# -4 \left|\begin{array}{}
# 2&4\\
# 3&2\end{array}\right|\\
# &=-5(18-10) +(12-20) - 4(4-12)\\
# &=-40 -8 + 32\\
# &=-16
# \end{align}

# * OK, so let's revisit the general determinant definition:
# \begin{equation}
#   D = |\arr{A}| = det(\arr{A})=\begin{cases}
#     a_{11} & \text{if $n=1$}.\\
#     a_{11}A_{11} + a_{12}A_{12} + ... + a_{1n}A_{1n} & \text{if n>1}.
#   \end{cases}
# \end{equation}
# where $A_{1k} = (-1)^{1+k} det(\arr{M_{1k}})$, for $k=1,2,...n$

# * So, for a 3 $\times$3 matrix, this would be 
# \begin{align}
# |\arr{A}| & = a_{11}(-1)^2 det(\arr{M_{11}}) + a_{12}(-1)^3 det(\arr{M_{12}}) + a_{13}(-1)^4 det(\arr{M_{13}})\\
# & =a_{11} det(\arr{M_{11}}) - a_{12} det(\arr{M_{12}}) + a_{13} det(\arr{M_{13}})
# \end{align}

# * This generalized definition is for when we work off of Row 1 only.
# * We can further generalize
# \begin{align}
# |\arr{A}| = a_{j1}A_{j1} + a_{j2}A_{j2} + ... + a_{jn}A_{jn} && \text{for n > 1}
# \end{align}
# where $A_{j1} = (-1)^{j+k}det(\arr{M_{jk}})$ (careful with the signs!)\
# $A_{j1}$ is the "cofactor" and $det(\arr{M_{jk}})$ is the "minor"

# * An n $\times$ n matrix has $n^2$ cofactors and $n^2$ minors.
# \begin{align}
# \begin{bmatrix}
# a_{11} & a_{12} & a_{13}\\
# a_{21} & a_{22} & a_{23}\\
# a_{31} & a_{32} & a_{33}
# \end{bmatrix} & \rightarrow 3 \times 3 & \text{will have 9 cofactors and 9 minors}\\
# \\
# \begin{array}{}
# A_{11}=+det(\arr{M_{11}}) && A_{12}=-det(\arr{M_{12}}) && A_{13}=+det(\arr{M_{13}})\\
# A_{21}=-det(\arr{M_{21}}) && A_{22}=+det(\arr{M_{22}}) && A_{23}=-det(\arr{M_{23}})\\
# A_{31}=+det(\arr{M_{31}}) && A_{32}=-det(\arr{M_{32}}) && A_{33}=+det(\arr{M_{33}})
# \end{array}
# \end{align}
# * Note the "checker board" of the signs:\
# \begin{array}{}
# + & - & +\\
# - & + & -\\
# + & - & +
# \end{array}
# * In general, working with Row 1 is easiest, but the super-generalized formula (allowing you to work with any row or column) can be useful depending on the matrix.\
# Ex: \begin{bmatrix}
# 9 & 2 & \bf{4}\\
# 5 & 3 & \bf{0}\\
# 1 & 6 & \bf{0}
# \end{bmatrix}\
# Which column or row to choose? Column 3, which has most zeros.
# \begin{align}
# |\arr{A}| &= a_{13}A_{13} + a_{23}A_{23} + a_{33}A_{33}\\
# & = a_{13} (-1)^{1+3} det(\arr{M_{13}})\\
# & = 4[(5)(6)-(3)(1)]\\
# & = 108
# \end{align}

# $\underline{Ques}:$ What if you have a row or column of all zeros?\
# $\underline{ANS}:$ Determinant is zero. We can "choose" all-zero row or column and then all $a_{jk}$'s will be zero.\
# Makes sense right? When $|\arr{A}|=0, \arr{A}$ is singular, or not invertible i.e. we can't convert the matrix to $\arr{I}$ if there is a row or column of all zeros.

# # In-class problem
# 
# Calculate the determinant of  
# \begin{align*}
# \begin{bmatrix}
# -1&2&3\\
# 0&5&3\\
# 3&1&1
# \end{bmatrix}
# \end{align*}
# 
# and upload it on gradescope under the in-class problem.

# 
# # Use of determinant in solving systems of equations
# 
# We will use the determinant a lot next class. Let's clarify a couple of uses now.
# 
# We're interested in a solution to the system $\arr{A}\vec{x}=\vec{b}$. There are two situations:
# * Non-homogeneous system of equations ($\vec{b}\neq\vec{0}$)
#   * If $det(\arr{A}\neq 0)$ or 
#   * rank($\arr{A}$)=n, or
#   * $\arr{A}$ is invertible,
#   * **THEN** a unique solution exists
# * Homogeneous system of equations ($\vec{b}=\vec{0}$)
#   * If $det(\arr{A}\neq 0)$
#     * There is only the trivial solution $\vec{x}=\vec{0}$
#   * If $det(\arr{A}= 0)$
#     * There is both the trivial solution $\vec{x}=\vec{0}$ **AND** a series of non-trivial solutions

# In[ ]:





# # Numerical calculation of determinant
# 
# Easy! np.linalg.det
# 
# Calculate the determinant of  
# \begin{align*}
# \begin{bmatrix}
# -1&2&3\\
# 0&5&3\\
# 3&1&1
# \end{bmatrix}
# \end{align*}
# 

# In[1]:


import numpy as np
A = np.array([[-1,2,3],[0,5,3],[3,1,1]])
np.linalg.det(A)


# In[ ]:




