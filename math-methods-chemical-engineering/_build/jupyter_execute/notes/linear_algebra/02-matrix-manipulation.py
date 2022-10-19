#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/ulissigroup/math-methods-chemical-engineering/blob/master/lecture_notes/02-matrix-manipulation.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$
# 

# # Note on HW expectations
# 
# This is an engineering course. The math and homeworks will difficult and challenging. I am expecting that they will take you on average about 8 hours per week to do (+ 4 hours in class, for a 12 credit class) including studying. If you find that they are taking you substantially longer to do, please let me know and we can discuss strategies. 

# # Systems of equations as matrices

# One use for a matrix is a set of linear equations!  $\arr{A}$: A set of linear equations
# 
# 
# Let's generalize this:
# \begin{align*}
# a_{11}x_1+a_{12}x_2+\dots+a_{1n}x_n&=b_1\\
# a_{21}x_1+a_{22}x_2+\dots+a_{2n}x_n&=b_2\\
# \vdots\\
# a_{m1}x_1+a_{m2}x_2+\dots+a_{mn}x_n&=b_m
# \end{align*}
# This is an $m \times n $ system of linear equations that can be represented by matrices and vectors!
# * $\arr{A}\vec{x}=\vec{b}$
# * We can also say  $\arr{A}=[a_{jk}]$, which just says that the matrix $\arr{A}$ consists of elements generalized as $a_{jk}$, where j=row number and k=column #
# * $[a_{jk}][x_k]=[b_j]$ - this is another notation that just says that we can multiple elements  of matrix $\arr{A}$ by the elements in $\vec{x}$ to get the elements of the vector $\vec{b}$, summing over the repeated index k.
# 
# 
# ## Examples
# 

# 
# ### Algebraic system
# We can take the system of two equations with three unknowns:
# \begin{align*}
# 2x_1+x_2+3x_3&=-1\\
# x_1-2x_2&=4
# \end{align*}
# and represent them as matrices and vectors:
# \begin{align*}
# \arr{A}\vec{x}&=\vec{b}\\
# \begin{bmatrix}
# 2&1&3\\
# 1&-2&0
# \end{bmatrix}
# \begin{bmatrix}
# x_1\\
# x_2\\
# x_3
# \end{bmatrix}&=
# \begin{bmatrix}
# -1\\4
# \end{bmatrix}
# \end{align*}
# 
# This uses matrix multiplication described in detail below!

# ### System of reactions
# 
# Let's consider a set of two reactions, and we know the rate of reaction of each:
# 
# 
# \begin{equation}
# \ce{2H2 + O2 ->[r_1] 2H2O}\\
# \ce{CO + H_2O ->[r_2] H2 + CO2}
# \end{equation}
# 
# Recall for your chemistry class that the stoichiometric coefficient is the integer in front of each species, and is positive for products and negative for reactants. 
# 
# We can represent this also as a system of equations!
# 
# \begin{align*}
# \arr{\alpha}\vec{r}&=\vec{N}\\
# \begin{bmatrix}
# -2&1\\
# -1&0\\
# 2&-1\\
# 0&-1\\
# 0&1
# \end{bmatrix}
# \begin{bmatrix}
# r_1\\
# r_2
# \end{bmatrix}&=
# \begin{bmatrix}
# r_{H2}\\
# r_{O2}\\
# r_{H2O}\\
# r_{CO}\\
# r_{CO2}
# \end{bmatrix}
# \end{align*}
# 
# We can use this to figure out how much of each species is being produced from the matrix of stoichiometric coefficients, and the vector of rates!

# # Fun ways to manipulate matrices

# ## Addition
# 
# 

# 
# If matrices/vectors are the same size, then you can add their elements together
# 
# ### Examples
# * $\arr{A}=
# \begin{bmatrix} 1&-1\\
# 0&-1\end{bmatrix}, \arr{B}=
# \begin{bmatrix} -2 & -1\\
# 2&4\end{bmatrix}$
#   * $\arr{A}+\arr{B}=\begin{bmatrix} -1 &-2\\2&3\end{bmatrix}$
# * $\vec{a}=\begin{bmatrix}5\\7\\-1\end{bmatrix}$, $\vec{b}=\begin{bmatrix}2\\-4\\0\end{bmatrix}$
#   * $\vec{a}+\vec{b}=\begin{bmatrix}7\\3\\-1\end{bmatrix}$
# 
# 
# 

# ### In-class Question
# 
# What is $\arr{A}+\vec{a}$?

#   * Undefined! 
#   * $\arr{A}$ is 2x2
#   * $\vec{a}$ is 3x1
#   * There is a size mismatch!

# 
# ## Transposition
# 
# 

# 
# If $\arr{A}=[a_{jk}]$, then $\arr{A}^T=[a_{kj}]$.
# 
# $^T$ is the transpose operation! In other words, we swap the columns with the rows.
# 
# 
# 
# #### Examples:
# * $\arr{A}=\begin{bmatrix}2&3\\1&4\end{bmatrix}$
#   * $\arr{A}^T=\begin{bmatrix}2&1\\3&4\end{bmatrix}$
# * $\vec{a}=\begin{bmatrix}-1&4&6\end{bmatrix}$ is a row vector (size 1x3)
#   * $\vec{a}^T=\begin{bmatrix}-1\\4\\6\end{bmatrix}$
# 
# We have to be really careful with the shapes when transposition things, since some operations (like addition) are only defined for arrays of certain sizes.
# 
# * $\vec{A}=\begin{bmatrix}2&4\\1&-2\\3&-5\end{bmatrix}$ is a **3x2 matrix**
#   * $\vec{A}^T=\begin{bmatrix}2&1&3\\4&-2&-5\end{bmatrix}$ is a **2x3 matrix**
# 
# ### Adding transposes
# 
# We aren't going to spend a lot of time proving various linear algebra identities, but some are really helpful. 
# 
# Consider $(\arr{A}+\arr{B})^T$
# * From the addition rules above, we know that the inner part is only defined if $\arr{A}$ and $\arr{B}$ are the same size.
# * Using the notation above $\arr{A}=[a_{jk}]$, $\arr{B}=[b_{jk}]$, $\arr{C}=\arr{A}+\arr{B}=[c_{jk}]$
#   * $(\arr{A}+\arr{B})^T=(a_{jk}+b_{jk})^T=(c_{jk})^T=c_{kj}=[a_{kj}]+[b_{kj}]=A^T+B^T$
# 
# #### Example:
# $\arr{A}=\begin{bmatrix} 1&0\\1&2\end{bmatrix},\arr{B}=\begin{bmatrix} -1&4\\0&5\end{bmatrix}$
# 
# \begin{align*}
# (\arr{A}+\arr{B})^T&=\arr{A}^T+\arr{B}^T\\
# \left(\begin{bmatrix} 1&0\\1&2\end{bmatrix}+\begin{bmatrix} -1&4\\0&5\end{bmatrix}\right)^T&=\begin{bmatrix} 1&0\\1&2\end{bmatrix}^T+\begin{bmatrix} -1&4\\0&5\end{bmatrix}^T\\
# \begin{bmatrix} 0&4\\1&7\end{bmatrix}^T&=\begin{bmatrix} 1&1\\0&2\end{bmatrix}+\begin{bmatrix} -1&0\\4&5\end{bmatrix}\\
# \begin{bmatrix} 0&1\\4&7\end{bmatrix}&=\begin{bmatrix} 0&1\\4&7\end{bmatrix}
# \end{align*}
# It worked. Amazing! Note that a single example of this working **is not a proof** of the identity!
# 
# 

# ## Matrix Multiplication
# 

# 
# $\arr{A}$: mxn matrix = $[a_{jk}]$
# 
# $\arr{B}$: rxp matrix = $[b_{jk}]$
# 
# The product $\arr{A}\arr{B}$ is defined **only if** n=r. That is, the number of columns in $\arr{A}$ must match the number of rows in $\arr{B}$.
# 
# **Definition:** $\arr{A}\arr{B}=\arr{C}$ is an mxp matrix, where
# \begin{align*}
# c_{jk}=\sum_{l=1}^na_{jl}b_{lk}
# \end{align*}
# 
# #### Examples: 
# 
# Say $\arr{A}=\begin{bmatrix}1&2&3\\ 0&-1&2\end{bmatrix}$, $\arr{B}=\begin{bmatrix}-2&1\\-1&1\\0&4\end{bmatrix}$
# 
# What is $\arr{C}=\arr{A}\arr{B}$?
# 

# * $\arr{A}$ is 2x3 and $\arr{B}$ is 3x2
#  * The product is defined because n=r=3
#  * the matrix $\arr{C}=\arr{A}\arr{B}$ is 2x2
#  * $c_{11}=a_{11}b_{11}+a_{12}b_{21}+a_{13}b_{31}\\
#      c_{12}=a_{11}b_{12}+a_{12}b_{22}+a_{13}b_{32}$ and so on
#  * So
#     *  $c_{11}=1\cdot(-2)+2\cdot(-1)+3\cdot 0=-4$
#     *  $c_{12}=1\cdot1+2\cdot1+3\cdot4=15$
#     *  $c_{21}=0\cdot(-2)+(-1)\cdot(-1)+2\cdot0=1$
#     *  $c_{22} = 0\cdot1 + (-1)\cdot 1+2\cdot 4=7$
#   * $\arr{C}=\begin{bmatrix}-4&15\\1&7\end{bmatrix}$
# 
# What is $\arr{C}=\arr{B}\arr{A}$?
# * B is a 3x2, A is a 2x3
# * The product is defined because 2=2
# * What is the size? 3x3!
# 

# 
# ### Properties of matrix multiplication
# * **Note:** This means that Matrix multiplication is not commutable! In general, $\arr{A}\arr{B}\neq \arr{B}\arr{A}$!!!!
# 
# * $\arr{A}\arr{B}=\arr{0}$ **does not** imply that $\arr{A}=\arr{0}, \arr{B}=0$, or that $\arr{B}\arr{A}=\arr{0}$
# * Matrix Multiplication is associative $\arr{A}(\arr{B}\arr{C})=(\arr{A}\arr{B})\arr{C}$, provided each product exists!
# * Matrix multiplication is distributive over addition:
#   * $(\arr{A}+\arr{B})\arr{C}=\arr{A}\arr{C}+\arr{B}\arr{C}\neq \arr{C}\arr{A}+\arr{C}\arr{B}=\arr{C}(\arr{A}+\arr{B})$

# # In-class exercise

# 
# For the following matrices:
# \begin{align*}
# \arr{D}=\begin{bmatrix}1&2&0\\0&3&-1\end{bmatrix}&&
# \arr{C}=\begin{bmatrix}1\\2\end{bmatrix}
# \end{align*}
# Calculate $\arr{D}^T\arr{C}$:
# * Predict the size of $\arr{D}^T\arr{C}$
# * Calculate $\arr{D}^T$
# * Calculate $\arr{D}^T\arr{C}$
# 
# Submit a picture of this on gradescope at the end of class!

# # Matrices in Python
# 
# * Many problems are too hard to solve by hand and must be solved numerically.
# * For these we use computational methods
# * We will extensively use Python to numerically solve problems in this course.
# * Why?
#   * Python is free
#   * You can use this anywhere you go
#   * Python does everything we need and much more
# * Python examples in these notes will be available to you through the syllabus
# * You should make sure you can run the examples, and that you get the same results
# * Ask questions when you do not understand
# 
# We'll do the python in-class example then come back to this

# ## Numpy arrays 
# 
# * Python doesn't know anything about arrays or linear algebra by default. Numpy "numerical python" is a library of standard objects and methods for numerical methods in python:
# https://docs.scipy.org/doc/numpy/user/quickstart.html . The quickstart is a little complicated but well-written.
# 
# * On top of basic numerical objects, scipy "scientific python" has additional routines for common operations we'll use, like integrating differential equations. More on that later!
# 
# * First, we have to import numpy to use

# In[1]:


import numpy as np

# Now we can use numpy. This makes an array of size 1 with the element 1
a = np.array([1])
print(a.shape)

#This is 1x1 2D array
a = np.array([[1]])
print(a.shape)

#Let's try a vector
a = np.array([1, 3, 5])
print(a.shape)


# In[2]:


a


# * NumPy’s main object is the homogeneous multidimensional array. It is a table of elements (usually numbers), all of the same type, indexed by a tuple of positive integers. In NumPy dimensions are called axes.
# 
# Let's try the in-class example.

# In[3]:


D = np.array([[1,2,0],
              [0,3,-1]])
C = np.array([[1,2]])

# This is a matrix multiply! D.T is the transpose of D
print(D.T@C)

#We can also do this with
print(np.matmul(D.T,C))

# The * is not matrix multiply. * is element-wise multiplication. Scarily, numpy
# will try to convert the second matrix to the right format by replicating it
# so it doesn't throw an error.

print(D.T*C)


# In[ ]:


D.shape


# In[ ]:


C = np.array([[1,2]])
C.shape


# In[ ]:


C = np.array([1,2])
C.shape


# Laika Count: 11

# In[ ]:




