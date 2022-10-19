#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/ulissigroup/math-methods-chemical-engineering/blob/master/lecture_notes/04-rank-gauss-examples.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
# 
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$
# 

# In[ ]:





# # Chemical Engineering example solved with Gauss elimination
# 
# Three tanks of water are attached in series. All tanks have the same cross-sectional area $A$. The flow rate through the valves is a function of the height of the water in the tanks (equations below)
# 
# <img src="https://docs.google.com/drawings/d/e/2PACX-1vRdeJhg87qBtUYmI28A_IHSRFfsJQUNfAiBTPTbcsiota93b2pL8MHzdpnXdhX21nG7nDAA1-mYLyh3/pub?w=632&amp;h=245">
# 
# 
# 

# ## Mass balances

# * Balance on tank 1 at steady state (accum=0)
# \begin{align*}
# \text{Accumulation}&=\text{In}-\text{Out}+\text{Generation}\\
# 0&=\rho F_0-\rho F_1\\
# 0&=F_0-F_1=F_0-\frac{h_1-h_2}{R_1}\\
# h_1-h_2&=R_1F_0
# \end{align*}
# To prepare for the augmented matrix, we'll keep unknowns on the left hand side, and constants on the right hand side!
# * Balance on tank 2 at steady state (accum=0)
# \begin{align*}
# \text{Accumulation}&=\text{In}-\text{Out}+\text{Generation}\\
# 0&=\rho\frac{h_1-h_2}{R_1}-\rho\frac{h_2-h_3}{R_2}=R_2h_1-R_2h_2-R_1h_2+R_1h_3\\
# R_2h_1-(R_1+R_2)h_2+R_1h_3&=0
# \end{align*}
# * Balance on tank 3 at steady state (accum=0)
# \begin{align*}
# \text{Accumulation}&=\text{In}-\text{Out}+\text{Generation}\\
# 0&=\rho\frac{h_2-h_3}{R_2}-\rho\frac{h_3}{R_3}=R_3h_2-R_3h_3-R_2h_3\\
# R_3h_2-(R_2+R_3)h_3&=0
# \end{align*}
# 
# 
# 

# ## Convert to math problem equation

# \begin{align*}
# \arr{A}\vec{h}=\vec{b}\\
# \begin{bmatrix}
# 1&-1&0\\
# R_2&-(R_1+R_2)&R_1\\
# 0&R_3&-(R_2+R_3)
# \end{bmatrix}
# \begin{bmatrix}
# h_1\\h_2\\h_3
# \end{bmatrix}&=\begin{bmatrix}
# R_1F_0\\0\\0
# \end{bmatrix}
# \end{align*}

# ## Solve system of equations

# For the situation where $A=5$ m$^2$, $F_0=5$ m$^3$/hr, $R_1=2$ hr/m$^2$, $R_2=1$ hr/m$^2$, $R_3=1$ hr/m$^2$, plug into matrix:
# \begin{align*}
# \begin{bmatrix}
# 1&-1&0\\
# 1&-3&2\\
# 0&1&-2
# \end{bmatrix}
# \begin{bmatrix}
# h_1\\h_2\\h_3
# \end{bmatrix}&=\begin{bmatrix}
# 10\\0\\0
# \end{bmatrix}
# \end{align*}
# 
# 
# 
# * Form $[\arr{A}|\vec{b}]$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&-1&0&10\\
# 1&-3&2&0\\
# 0&1&-2&0
# \end{array}\right]
# \end{align*}
# * $R_2=R_2-R_1$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&-1&0&10\\
# 0&-2&2&-10\\
# 0&1&-2&0
# \end{array}\right]
# \end{align*}
# * $R_2=-R_2/2$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&-1&0&10\\
# 0&1&-1&5\\
# 0&1&-2&0
# \end{array}\right]
# \end{align*}
# * $R_3=R_2-R_3$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&-1&0&10\\
# 0&1&-1&5\\
# 0&0&1&5
# \end{array}\right]
# \end{align*}
# * $R_2=R_2+R_3$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&-1&0&10\\
# 0&1&0&10\\
# 0&0&1&5
# \end{array}\right]
# \end{align*}
# * $R_1=R_1+R_2$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&0&0&20\\
# 0&1&0&10\\
# 0&0&1&5
# \end{array}\right]
# \end{align*}
# 
# Math solved!

# ## Solve engineering problem

# * $h_1=20$. What units? Remember $h\sim FR=(m^3/hr)(hr/m^2)=m$
# * $h_1=20$ m, $h_2=10$ m, $h_3=5$ m.
# 
# Can this solution make physical sense? Yes!

# # Rank of a matrix

# * The rank of a matrix is the # of non-zero rows after gauss elimination
# * This tells us about the linear independent of the our system
# 

# ## Example

# \begin{align*}
# \arr{A}=\begin{bmatrix}
# -3&1&-1\\
# 1&0&1\\
# -2&2&2
# \end{bmatrix}
# \end{align*}
# * After Gauss elimination
# \begin{align*}
# \arr{A}=\begin{bmatrix}
# 1&0&1\\
# 0&1&2\\
# 0&0&0
# \end{bmatrix}
# \end{align*}
# *There are two nonzero rows, so rank$(\arr{A})=2$. 
# * This also means that only two of these equations are linearly independent!
# * The reduced echelon matrix tells us that $x_1+x_3$ and $x_1+2x_3$ form **a basis** of the system!
# 

# ## Example 2

# Let's say we get an augment matrix after gauss elimination of:
# \begin{align*}
# \arr{\tilde{A}}=\begin{bmatrix}
# -3&1&-1&-3\\
# 1&0&1&5\\
# -2&2&2&-17
# \end{bmatrix}
# \end{align*}
# * rank$(\arr{\tilde{A}})=3$
# * rank$(\arr{A})=2$
# 
# **No solution!!**
# 
# For an m=n system like the 3x3 system here:
# * If rank$(\arr{A})=$rank$(\arr{\tilde{A}})=n$, one unique solution!
# * If rank$(\arr{A})<$rank$(\arr{\tilde{A}})=n$, no solution!

# # Inner (dot) products of vectors
# 
# You've probably done dot products in your calculus class, but let's go ahead and review dot products.
# 
# * This is a special case of matrix multiplication that occurs very often
# * If $\vec{a}$ and $\vec{b}$ are vectors, both with the same number of elements, then matrix multiplication yields a scalar (a 1x1 matrix) called the inner or dot  product
# \begin{align*}
# \vec{a}\cdot\vec{b}=\vec{a}^T\vec{b}=\begin{bmatrix}a_1&\dots&a_n\end{bmatrix}
# \begin{bmatrix}b_1\\\vdots\\b_n\end{bmatrix}=\sum_{l=1}^n a_lb_l
# \end{align*}
# 
# ## Example
# 
# \begin{align*}
# \vec{a}=\begin{bmatrix}
# 4\\-1\\2
# \end{bmatrix}, \vec{b}=\begin{bmatrix}
# 1\\0\\3
# \end{bmatrix}\\
# \vec{a}\cdot\vec{b}=4(1)+-1(0)+2(3)=10
# \end{align*}
# 
# * Matrix multiplication is just a bunch of dot products. For example, in $\arr{C}=\arr{A}\arr{B}$, every operation is just a dot product of a row and columnn vector

# # Diagonal Matrices

# * Definition: Square matrices that have non-zero entries only on the diagonal. Any entry above of below the diagonal is zero
# * Example:
# \begin{align*}
# \begin{bmatrix}
# 0.2&0&0\\0&5&0\\0&0&1
# \end{bmatrix}
# \end{align*}

# ## Identity (unit) matrix
# 
# * A square matrix  $\arr{I}$ where every element in the diagonal is 1, and all other elements are 0. 
# \begin{align*}
# \arr{I}=\begin{bmatrix}1&0\\0&1\end{bmatrix} &&\arr{I}=\begin{bmatrix}1&0&0\\0&1&0\\0&0&1\end{bmatrix} &&etc
# \end{align*}
# * We can denote the size with $\arr{I}_n$
# * **Special property!** Multiplying a matrix $\arr{A}$ by the same size identity matrix $\arr{I}$ yields back $\arr{A}$!
# \begin{align*}
# \arr{I}\arr{A}=\arr{A}\arr{I}=\begin{bmatrix}
# a_{11}&a_{12}\\
# a_{21}&a_{22}
# \end{bmatrix}\begin{bmatrix}
# 1&0\\0&1
# \end{bmatrix}=\begin{bmatrix}
# a_{11}&a_{12}\\
# a_{21}&a_{22}
# \end{bmatrix}
# \end{align*}

# ## Scalar matrix
# 
# * All non-zero entries of a diagonal matrix $\arr{S}$ are equal to  some number $c$. In other words, $\arr{S}=c\arr{I}$
# \begin{align*}
# \arr{S}=\begin{bmatrix} c&0\\0&c\end{bmatrix} &&e.g. && \arr{S}=\begin{bmatrix} 3&0\\0&3\end{bmatrix}
# \end{align*}
# * From the unit identity info above, $\arr{S}\arr{A}=c\arr{I}\arr{A}=\arr{A}\arr{S}=\arr{A}c\arr{I}=c\arr{A}$

# # Matrix inversion
# 
# * A real number $a$ has a multiplicative inverse if there exists a $b$ such that $ab=1$
# * Any non-zero $a$ has a multiplicative inverse $b=1/a$.
# 
# We can generalize this idea to matrices!
# 
# * A square matrix $\arr{A}$ is **invertible** or **non-singular** if there is a matrix $\arr{B}$ such that
# \begin{align*}
# \arr{A}\arr{B}=\arr{B}\arr{A}=\arr{I}
# \end{align*}
#   * Ex: $\arr{A}=\begin{bmatrix}2&4\\3&1\end{bmatrix}$ and $\arr{B}=\begin{bmatrix}-1/10&2/5\\3/10&-1/5\end{bmatrix}$ are inverses of each other because $\arr{A}\arr{B}=\arr{B}\arr{A}=\arr{I}$
# * We say $\arr{B}$ is the multiplicative inverse of $\arr{A}$
# * If $\arr{B}$ and $\arr{C}$ are both inverses of $\arr{A}$, then 
# \begin{align*}
# \arr{B}=\arr{B}\arr{I}=\arr{B}(\arr{A}\arr{C})=(\arr{B}\arr{A})\arr{C}=\arr{I}\arr{C}=\arr{C} \\\therefore \arr{B}=\arr{C} \text{ and $\arr{A}$ has at most one inverse}
# \end{align*}
# * The inverse of $\arr{A}$ is denoted $\arr{A}^{-1}$
# * A square matrix is called **singular** if it doesn't have an inverse
# 
# 
# 
# 

# ## Computing an inverse

# We already have the tools to compute the inverse!! 
# 
# * For a square matrix $\arr{A}$, we will augment it with the identity matrix
# * We we will use Gauss elimination to convert $\arr{A}$ to the identity matrix
# * The augmented matrix is then the inverse! 
# * Magic! (not really)

# ### Example

# Find $\arr{A}^{-1}$ for $\arr{A}=\begin{bmatrix}2&-6\\4&-2\end{bmatrix}$
# 
# * Form augmented matrix
# \begin{align*}
# \left[\begin{array}{rr|rr}
# 2&-6&1&0\\
# 4&-2&0&1
# \end{array}\right]
# \end{align*}
# * $R_1=R_1/2$
# \begin{align*}
# \left[\begin{array}{rr|rr}
# 1&-3&1/2&0\\
# 4&-2&0&1
# \end{array}\right]
# \end{align*}
# * $R_2=R_2-4R_1$
# \begin{align*}
# \left[\begin{array}{rr|rr}
# 1&-3&1/2&0\\
# 0&10&-2&1
# \end{array}\right]
# \end{align*}
# * $R_2=R_2/10$
# \begin{align*}
# \left[\begin{array}{rr|rr}
# 1&-3&1/2&0\\
# 0&1&-1/5&1/10
# \end{array}\right]
# \end{align*}
# * $R_1=R_1+3R_2$
# \begin{align*}
# \left[\begin{array}{rr|rr}
# 1&0&-1/10&3/10\\
# 0&1&-1/5&1/10
# \end{array}\right]
# \end{align*}
# 
# * $\arr{A}^{-1}=\begin{bmatrix}-1/10&3/10\\
# -1/5&1/10
# \end{bmatrix}$
# * Quick check!
# \begin{align*}
# \arr{A}\arr{A}^{-1}=\arr{I}=\arr{A}^{-1}\arr{A}
# \end{align*}
# 
# Let's do this one with python!
# 

# In[1]:


import numpy as np

#Define A and the inverse of A
A = np.array([[2,-6],[4,-2]])
Ainv = np.array([[-1/10,3/10],[-1/5,1/10]])

#Note that this looks very good, but some residual rounding errors in the 
# computation has one element VERY close to zero but not quite right. For 
# practical purposes, anything less than about 10^-8 is probably the same 
# thing as zero in numerical methods
print(A@Ainv)

#Check Ainv*A
print(Ainv@A)


# # Rank, inverses, and row echelon form, and solving systems of linear equations in numpy

# ## Rank of a matrix

# The ndim property of a numpy array is the same thing as the rank. For the example above:

# In[2]:


A = np.array([[-3,1,-1],[1,0,1],[-2,2,2]])
print(A)

# We can use np.linalg.matrix_rank(A)
print('The rank of A is %d'%np.linalg.matrix_rank(A))


# Let's adjust the matrix and see if the rank changes ($a_{22}=2$ to $a_{22}=3$)

# In[3]:


A = np.array([[-3,1,-1],[1,0,1],[-2,2,3]])

# Calculate the matrix rank of A
print('The rank of A is %d'%np.linalg.matrix_rank(A))


# Confusingly, np.rank(A) **is not** the matrix rank of A.

# ## Inverse of a matrix
# 
# Inverses are also easy with np.linalg.inverse

# In[4]:


A = np.array([[2,-6],[4,-2]])

#Same as the answer we got above!
print(np.linalg.inv(A))


# We can ask numpy what the inverse of a singular matrix is, and it will give us an answer, but that answer is garbage.

# In[5]:


#Try to get the inverse of a singular matrix
A = np.array([[-3,1,-1],[1,0,1],[-2,2,2]])

#Same as the answer we got above!
Ainv = np.linalg.inv(A) 
print(Ainv)

#Check; notice that Ainv*A is not the identity! This is a problem with numerical
# methods; we have to be on top of our game all the time!
Ainv@A


# In[6]:


Ainv


# ## Gauss elimination
# 
# * There's not a good method to do Gauss elimination with numerical methods. You can do it you use symbolic math packages (sort of like mathematica, but in python), but these largely defeat the point of using numerical methods! 
# 
# * The closest thing you can do is get a matrix into upper diagonal form, which is echelon form, but this only works for square matrics and isn't really helpful (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.lu_factor.html)
# 
# * We use Gauss elimination to solve linear systems and calculate inverses, and there are better ways to do these. 

# ## Solving a system of linear equations
# 
# * There are excellent methods to solve a system of linear equations (https://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.solve.html) that is determined (square) and full rank
# * solve takes in the matrix $\arr{A}$ and the vector $\vec{b}$ in $\arr{A}\vec{x}=\vec{b}$ and returns $\vec{x}$
# 

# ### Example for determined system
# 
# 
# Let's try the example from above for the three-tank system
# 
# \begin{align*}
# \begin{bmatrix}
# 1&-1&0\\
# 1&-3&2\\
# 0&1&-2
# \end{bmatrix}
# \begin{bmatrix}
# h_1\\h_2\\h_3
# \end{bmatrix}&=\begin{bmatrix}
# 10\\0\\0
# \end{bmatrix}
# \end{align*}

# In[7]:


A = np.array([[1,-1,0],
              [1,-3,2],
              [0,1,-2]])
b = np.array([10,0,0])

#Same answer as above!
print(np.linalg.solve(A,b))


# This works really well for a determined system.

# ### Example for an underdetermined system
# 
# Let's do one of the examples from the previous class
# 
# \begin{align*}
# x_1+2x_2+x_3=1\\
# 2x_1-x_2+x_2=2\\
# 4x_1+3x_2+3x_3=4\\
# 3x_1+x_2+2x_3=3
# \end{align*}

# In[8]:


A = np.array([[1,2,1],
              [2,-1,1],
              [4,3,3],
              [3,1,2]])
b = np.array([1,2,4,3])


#Same answer!
print(np.linalg.solve(A,b))


# Notice that np.linalg.solve can't handle this case. We can still get a solution,with np.linalg.lstsq, which tries to find a vector that satisfies this. 

# In[ ]:


#Try least squares solution to this problem
x, res, rank, s = np.linalg.lstsq(A,b)

#First, print the rank of the matrix:
print(rank)

#print the solution
print(x)

#Check that this is a solution
print(A@x)


# So, np.linalg.lstsq does give us **a solution** but not **all of the solutions**. There are other solutions that also would have worked as we saw. The simplest is [1,0,0] 

# In[ ]:


A@np.array([1,0,0])


# Once again, numerical methods can work well, but we have to be extra sharp to make sure we understand what they are doing!

# # In-class assignment
# 
# For the example we did at the beginning of class today (the three tank problem), set up the matrix A and vector b, check the rank of the matrix A, and then solve for the solution x. 

# Laika count 40; check the np.linalg.matrix_rank example

# In[ ]:




