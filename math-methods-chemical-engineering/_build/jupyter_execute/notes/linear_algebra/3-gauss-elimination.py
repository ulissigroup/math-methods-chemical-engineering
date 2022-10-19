#!/usr/bin/env python
# coding: utf-8

# {{ badge }}
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$
# 

# # Return to Systems of Equations

# Now that we understand matrix multiplication, let's go back to linear systems of equations:
# \begin{align*}
# \begin{bmatrix}
# a_{11}&a_{12}&\dots&a_{1n}\\
# a_{21}&a_{22}&\dots&a_{2n}\\
# \vdots&\vdots&\vdots&\vdots\\
# a_{m1}&a_{m2}&\dots&a_{mn}
# \end{bmatrix}
# \begin{bmatrix}
# x_1\\x_2\\\vdots\\x_n
# \end{bmatrix}=
# \begin{bmatrix}
# b_1\\b_2\\\vdots\\b_m
# \end{bmatrix}
# \end{align*}
# 
# This is an $mxn$ matrix multiplied by a nx1 vector, resulting in an mx1 vector. An example is the same system of equations we've solved before:
# \begin{align*}
# \begin{bmatrix}
# 2&1&3\\
# 1&-2&0
# \end{bmatrix}
# \begin{bmatrix}
# x_1\\x_2\\x_3
# \end{bmatrix}=
# \begin{bmatrix}
# -1\\4
# \end{bmatrix}\\
# \begin{bmatrix}
# 2x_1+x_2+3x_3\\
# x_1-2x_2+0x_3
# \end{bmatrix}=
# \begin{bmatrix}
# -1\\4
# \end{bmatrix}
# \end{align*}
# 
# 
# We've talked a lot about how we were going to use linear algebra to solve systems of equations, so it's time to actually solve them!

# # Gauss Elimination
# 

# 
# This is a standard method for solving linear systems. It is not the only way you can do it. Let's do an example without matrices first:
# \begin{align*}
# 2x_1+5x_2=2\\
# 4x_1+3x_2=18
# \end{align*}
# or 
# \begin{align*}
# \begin{bmatrix}
# 2&5\\4&3
# \end{bmatrix}\begin{bmatrix}
# x_1\\x_2
# \end{bmatrix}&=
# \begin{bmatrix}
# 2\\18
# \end{bmatrix}
# \end{align*}
# 
# * Let's first solve for one variable, then solve for the other using back substitution
# * One approach is to multiple first equation by 2, then subtract the second equation:
# \begin{align*}
# 4x_1+10x_2=4\\
# -[4x_1+3x_2=18]\\
# \rule{4cm}{0.4pt}\\
# \rightarrow 7x_2=-14\\
# x_2=-2
# \end{align*}
# * We can now use this result to find that $x_1=6$. After some operations, we arrived at the system of equations:
# \begin{align*}
# x_1&=6\\
# x_2&=-2\\
# \begin{bmatrix}
# 1&0\\0&1
# \end{bmatrix}\begin{bmatrix}
# x_1\\x_2
# \end{bmatrix}&=
# \begin{bmatrix}
# 6\\-2
# \end{bmatrix}
# \end{align*}
# * This approach involved manipulation of whole equations, one at a time.

# ## Elementary Operation for Equations
# 
# These are allowed when solving equations:
# 
# 1.   Interchange of two equations (no effect!)
# 2.   Addition or subtraction of one equation with another
# 3. Multiplication of an equation by a non-zero constant. Q: Why is this only for non-zero?
# 
# This is essentially the rules for Gauss elimination but with equations instead of matricies. Let's do them side by side. The same elementary operations apply to matrices, just substitute "equation" for "row"!
# 

# ### Examples
# 
# Standard system of linear system of equations: 
# \begin{align*}
# \begin{array}{r}
# x_1+2x_2-x_3=4\\
# 4x_2-2x_3=-2\\
# x_1-2x_2+3x_3=0
# \end{array}
# \end{align*}
# 
# 
# 

# #### Form the augmented matrix
# 
# First, we represent this as $\arr{A}$ and $\vec{b}$
# \begin{align*}
# \arr{A}=\begin{bmatrix} 
# 1&2&-1\\
# 0&4&-2\\
# 1&-2&3
# \end{bmatrix} && \vec{b}= 
# \begin{bmatrix}
# 4\\
# -2\\
# 0
# \end{bmatrix}
# \end{align*}
# 
# We then form an augmented matrix $\arr{\tilde{A}}=[\arr{A}:\vec{b}]$:
# \begin{align*}
# \begin{array}{r}
# x_1+2x_2-x_3=4\\
# 4x_2-2x_3=-2\\
# x_1-2x_2+3x_3=0
# \end{array}&& \left[\begin{array}{rrr|r}
# 1&2&-1&4\\
# 0&4&-2&-2\\
# 1&-2&3&0
# \end{array}\right]
# \end{align*}
# 

# 
# #### Augmented matrix to upper triangular form
# 
# 

# 
# Example procedure for Gauss Elimination 
# * Let's swap rows two and three
# \begin{align*}
# \begin{array}{r}
# x_1+2x_2-x_3=4\\
# x_1-2x_2+3x_3=0\\
# 4x_2-2x_3=-2
# \end{array}&& \left[\begin{array}{rrr|r}
# 1&2&-1&4\\
# 1&-2&3&0\\
# 0&4&-2&-2
# \end{array}\right]
# \end{align*}
# 
# * Subtract row 1 from row 2 and replace row 2
# \begin{align*}
# \begin{array}{r}
# x_1+2x_2-x_3=4\\
# -4x_2+4x_3=-4\\
# 4x_2-2x_3=-2
# \end{array}&& \left[\begin{array}{rrr|r}
# 1&2&-1&4\\
# 0&-4&4&-4\\
# 0&4&-2&-2
# \end{array}\right]
# \end{align*}
# 
# *  Add rows 2 and 3 to replace row 3
# \begin{align*}
# \begin{array}{r}
# x_1+2x_2-x_3=4\\
# -4x_2+4x_3=-4\\
# 2x_3=-6
# \end{array}&& \left[\begin{array}{rrr|r}
# 1&2&-1&4\\
# 0&-4&4&-4\\
# 0&0&2&-6\end{array}\right]
# \end{align*}
# 
# 
# **This matrix is now in upper triangular form! Only the diagonal and things above it are occupied. Notice it would now be pretty easy to solve this system.**

# #### Upper triangular to row echelon form

# *  Divide row 3 by 2 and row 2 by -4
# \begin{align*}
# \begin{array}{r}
# x_1+2x_2-x_3=4\\
# x_2-x_3=1\\
# x_3=-3
# \end{array}&& \left[\begin{array}{rrr|r}
# 1&2&-1&4\\
# 0&1&-1&1\\
# 0&0&1&-3\end{array}\right]
# \end{align*}
# 
# **The matrix is now in "row echelon" or "echelon" form!**
# * The first non-zero entry of each row is 1
# * Any rows of all zeros must be at the bottom of the matrix
# * Rows with more leading zeros must be underneath those with less
# 
# We could solve this pretty easily but it  would still take some algebra. Almost there! 

# #### Echelon to reduced row echelon form

# *  Add rows 2&3 to replace row 2
# \begin{align*}
# \begin{array}{r}
# x_1+2x_2-x_3=4\\
# x_2=-2\\
# x_3=-3
# \end{array}&& \left[\begin{array}{rrr|r}
# 1&2&-1&4\\
# 0&1&0&-2\\
# 0&0&1&-3\end{array}\right]
# \end{align*}
# *  Add rows 1&3 to replace row 1
# \begin{align*}
# \begin{array}{r}`
# x_1+2x_2=1\\
# x_2=-2\\
# x_3=-3
# \end{array}&& \left[\begin{array}{rrr|r}
# 1&2&0&1\\
# 0&1&0&-2\\
# 0&0&1&-3\end{array}\right]
# \end{align*}
# *  Subtract 2x row 2 from row 1 to replace row 1
# \begin{align*}
# \begin{array}{r}
# x_1=5\\
# x_2=-2\\
# x_3=-3
# \end{array}&& \left[\begin{array}{rrr|r}
# 1&0&0&5\\
# 0&1&0&-2\\
# 0&0&1&-3\end{array}\right]
# \end{align*}
# 
# The "solved" matrix is now in "reduced row echelon form" or "reduced echelon" form. 
# * Satisfies the definition of echelon form
# * **AND** the 1st non-zero entry in each row is the only non-zero entry in each column!
# 

# ## Summary of steps

# 1. Convert equations into an augmented matrix
# 2. Put the matrix in upper triangular form
#   * Start by converting element $a_{11}$ to a 1
#   * Use row 1 to convert elements $a_{j1}$ to zeros
#   * Convert any other elements below the diagonal to zero if possible
# 3. Put matrix in echelon form by converting diagonal elements to 1
# 4. Put matrix in reduced row echelon form by converting as many elements above the diagonal as possible to zero 
# 

# ## In-class example
# 
# Convert the following system of equations to reduced row echelon form
# 
# \begin{align*}
# 2x_1+4x_2&=-4\\
# 3x_1-x_2&=8
# \end{align*}
# 

# ### Solution
# 
# * Form $[\arr{A}|\vec{b}]$
# \begin{align*}
# \left[\begin{array}{rr|r}
# 2&4&-4\\
# 3&-1&8
# \end{array}\right]
# \end{align*}
# * $R_1=1/2R_1$
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&2&-2\\
# 3&-1&8
# \end{array}\right]
# \end{align*}
# * $R_2=R_2-3R_1$ **Upper triangular!**
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&2&-2\\
# 0&-7&14
# \end{array}\right]
# \end{align*}
# * $R_2=-1/7R_2$ **Echelon form!**
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&2&-2\\
# 0&1&-2
# \end{array}\right]
# \end{align*}
# * $R_1=R_1-2R_2$ **Reduced echelon!**
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&0&2\\
# 0&1&-2
# \end{array}\right]
# \end{align*}
# 
# 
# Final solution:
# 
# \begin{align*}
# x_1=2\\
# x_2=-2
# \end{align*}
# 
# Last thing: check that this solution works for each equation!

# ## Gauss elimination as a test for linear independence
# 

# 
# In addition to solving sets of linear equations, Gauss elimination is a powerful way to look for **linear independence**
# * For example, $x_1+3x_2=-1$ and $2x_1+6x_2=-2$ are not linearly independent. That is, they both contain the same information!
# 
# * For equations and variables, systems can be:
#   * **Underdetermined:** Less equations than unknowns (m<n)
#   * **Overdetermined:** More equations that unknowns (m>n)
#   * **Determined** Same # of equations as unknowns (m=n)
# 
# * Gauss elimination will work for all of these scenarios
# * Possible solution scenarios are:
#   * Infinite solutions
#   * A unique solution
#   * No solution
# 
# * For two unknowns, we can draw this as lines:
# 
# <img src="https://docs.google.com/drawings/d/e/2PACX-1vSfDcBCxRtHnHRryC0HlWY51Tsad1g8QmZ0BQj-ksCBpGUtNjGIL52eyNa52huXKzL7T0PWHtHKe1UL/pub?w=659&amp;h=249">
# 
# * For two unknowns, we either have 0, 1, or infinited solutions (can't have 2 solutions for example)
# 
# 

# ### Examples for underdetermined systems

# #### Example 1
# 
# \begin{align*}
# x_1+x_2+x_3+x_4+x_5&=2\\
# x_1+x_2+x_3+2x_4+2x_5&=3\\
# x_1+x_2+x_3+2x_4+3x_5&=2
# \end{align*}
# 
# Three equations, five unknowns. Suspicious!

#  Let's get this to reduced echelon form:
# 
# 
# * Form $[\arr{A}|\vec{b}]$
# \begin{align*}
# \left[\begin{array}{rrrrr|r}
# 1&1&1&1&1&2\\
# 1&1&1&2&2&3\\
# 1&1&1&2&3&2
# \end{array}\right]
# \end{align*}
# * $R_2=R_2-R_1$, $R_3=R_3-R_1$
# \begin{align*}
# \left[\begin{array}{rrrrr|r}
# 1&1&1&1&1&2\\
# 0&0&0&1&1&1\\
# 0&0&0&1&2&0
# \end{array}\right]
# \end{align*}
# * $R_3=R_3-R_2$ **Echelon form!**
# \begin{align*}
# \left[\begin{array}{rrrrr|r}
# 1&1&1&1&1&2\\
# 0&0&0&1&1&1\\
# 0&0&0&0&1&-1
# \end{array}\right]
# \end{align*}
# Notice that there are extra blocks of zeros in row 2/3. This will lead to infinite solutions!
# * $R_1=R_1-R_3$, $R_2=R_2-R_3$
# \begin{align*}
# \left[\begin{array}{rrrrr|r}
# 1&1&1&1&0&3\\
# 0&0&0&1&0&2\\
# 0&0&0&0&1&-1
# \end{array}\right]
# \end{align*}
# * $R_1=R_1-R_2$ **Reduced echelon!**
# \begin{align*}
# \left[\begin{array}{rrrrr|r}
# 1&1&1&0&0&1\\
# 0&0&0&1&0&2\\
# 0&0&0&0&1&-1
# \end{array}\right]
# \end{align*}
# * This corresponds to the system of equations
# \begin{align*}
# x_1+x_2+x_3=1\\
# x_4=2\\
# x_5=-1
# \end{align*}
# * This means that the solution vector 
# \begin{align*}
# \vec{x}=\begin{bmatrix}
# \alpha\\ \beta\\ 1-\alpha-\beta\\ 2\\-1
# \end{bmatrix}
# \end{align*}
# 
# **for any values $\alpha, \beta$!!**. 
# * There are an infinite number of solutions, which is usually but not always the case with an underdetermined system!
# 

# #### Example 2
# 
# Let's do another example!
# 
# \begin{align*}
# x_1+2x_2+x_3&=1\\
# 2x_1+4x_2+2x_3&=3
# \end{align*}

# 
# 
# * Two equations, three unknowns
# * Form $[\arr{A}|\vec{b}]$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&2&1&1\\
# 2&4&2&3
# \end{array}\right]
# \end{align*}
# * $R_2=R_2-2R_1$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&2&1&1\\
# 0&0&0&1
# \end{array}\right]
# \end{align*}
# * $0\neq 1$, so there is no solution and this is an inconsistent system!!!

# ### Overdetermined systems

# #### Example 1
# 
# 
# \begin{align*}
# x_1+2x_2=1\\
# x_1+4x_2=5\\
# x_1+3x_2=3
# \end{align*}
# 

# 
# * Three equations, two unknowns. Suspicious!
# * Form $[\arr{A}|\vec{b}]$
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&2&1\\
# 1&4&5\\
# 1&3&3
# \end{array}\right]
# \end{align*}
# * $R_2=R_2-R_1$, $R_3=R_3-R_1$
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&2&1\\
# 0&2&4\\
# 0&1&2
# \end{array}\right]
# \end{align*}
# * $R_3=R_3-1/2R_2$, $R_2=1/2R_2$
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&2&1\\
# 0&1&2\\
# 0&0&0
# \end{array}\right]
# \end{align*}
# Note the last equation doesn't really have any information! 
# * $R_1=R_1-2R_2$
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&0&-3\\
# 0&1&2\\
# 0&0&0
# \end{array}\right]
# \end{align*}
# 
# * The solution is thus $x_1=-3, x_2=2$. One unique solution!

# #### Example 2
# 
# \begin{align*}
# x_1+2x_2+x_3=1\\
# 2x_1-x_2+x_3=2\\
# 4x_1+3x_2+3x_3=4\\
# 3x_1+x_2+2x_3=3
# \end{align*}

# 
# * Four equations, three unknowns. Suspicious!
# * Form $[\arr{A}|\vec{b}]$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&2&1&1\\
# 2&-1&1&2\\
# 4&3&3&4\\
# 3&1&2&3
# \end{array}\right]
# \end{align*}
# * $R_2=R_2-2R_1$, $R_3=R_3-4R_1$, $R_4=R_4-3R_1$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&2&1&1\\
# 0&-5&-1&0\\
# 0&-5&-1&0\\
# 0&-5&-1&0
# \end{array}\right]
# \end{align*}
# * $R_2=-1/5R_2$, $R_3=R_3-R_2$, $R_4=R_4-R_2$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&2&1&1\\
# 0&1&1/5&0\\
# 0&0&0&0\\
# 0&0&0&0
# \end{array}\right]
# \end{align*}
# * $R_1=R_1-2R_2$
# \begin{align*}
# \left[\begin{array}{rrr|r}
# 1&0&3/5&1\\
# 0&1&1/5&0\\
# 0&0&0&0\\
# 0&0&0&0
# \end{array}\right]
# \end{align*}
# 
# This corresponds to
# \begin{align*}
# x_1+\frac{3}{5}x_3&=1\\
# x_2+\frac{1}{5}x_3&=0
# \end{align*} 
# The solution can be written as
# \begin{align*}
# \vec{x}=\begin{bmatrix}
# 1-\frac{3}{5}\alpha\\
# \frac{-1}{5}\alpha\\
# \alpha
# \end{bmatrix}
# \end{align*}
# 
# **Infinite number of solutions!!**

# #### Example 3
# 
# 
# \begin{align*}
# x_1+x_2=1\\
# x_1-x_2=3\\
# -x_1+2x_2=-2
# \end{align*}
# 

# 
# * Three equations, two unknowns. Suspicious!
# * Form $[\arr{A}|\vec{b}]$
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&1&1\\
# 1&-1&3\\
# -1&2&-2
# \end{array}\right]
# \end{align*}
# * $R_2=R_2-R_1$, $R_3=R_3+R_1$
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&1&1\\
# 0&-2&2\\
# 0&3&-1
# \end{array}\right]
# \end{align*}
# * $R_2=\frac{-1}{2}R_2$, $R_3=R_3/3$
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&1&0\\
# 0&1&-1\\
# 0&1&-1/3
# \end{array}\right]
# \end{align*}
# * $R_3=R_3-R_2$
# \begin{align*}
# \left[\begin{array}{rr|r}
# 1&1&1\\
# 0&1&-1\\
# 0&0&-2/3
# \end{array}\right]
# \end{align*}
# * $0\neq -2/3$. **Inconsistent, no solutions!!**

# ## Final notes
# 
# * Underdetermined systems
#   * Often infinite solutions
#   * Can be infinite, no solution, or one solution
# * Over-determined system
#   * Often no solution
#   * Can be infinite, no solution, or one solution
# * Determined system
#   * Often one solution
#   * Can be infinite, no solution, or one solution
# 
# **Homogeneous systems (when $\vec{b}=0$) always has at least one solution! ($\vec{x}=\vec{0}$)**

# In[ ]:




