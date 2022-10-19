#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/ulissigroup/math-methods-chemical-engineering/blob/master/lecture_notes/01-intro-matrices.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$    
# $$\require{mhchem}$$
# 

# # Introduction
# 
# * The purpose of this course is to learn how to formulate and solve chemical engineering problems using analytical and numerical math methods.
# 
# * This is both a math course and an engineering course. I will teach you mathematical tools to solve real world problems.
# 
# * The best chemical engineers are often the ones that are best at math. Math is a powerful tool and chemical engineers are jacks of many trades that requires a hefty toolbox.
# 
# * Being good at math is not some talent that you were born with. Every one of you has the opportunity to be good at math this semester and for the rest of your life. All you need to do is practice. A lot.
# 
# * There will be many demands on your time this semester. Lab 1, Fluids, electives, extracurriculars. We ask a lot of you because we know you can do it! Working hard in this class is an investment that will help you with the rest of your classes. 

# # How to solve a generic engineering problem
# 
# 1. **Problem Statement**
#   * Design a reactor; design a vehicle for drug delivery
# 2.  **Problem formulation** 
#   * convert the statement into a math model if possible (i.e. ill-posed problems exist)
# 3. **Problem analysis**
#   * Identify the type of math necessary to solve the problem, and determine a reasonable approach (analytical or numerical).
# 4. **Problem solution**
#   * solve the **mathematical** problem (which is not the same as the engineering problem)
# 5. **Problem evaluation**
#   * Sanity check. Does the math solution make physical sense? Does it help you solve the engineering problem? If not, 1/2/3/4 may need revision.
# 
# Most of your math courses only dealt with step #3. But Chemical Engineers need to be good at all of the steps. That's why we have this class to teach you the math tools and practice all of the steps together. 
# 
# 
# 

# # Methods in this course
# 
# This course focused on 4 main methods:
# 1. Linear Algebra - e.g. solving systems of linear equations
# 2. Differential Equations - govern rates of change (space or time) of physical quantitities
# 3. Probability & Stats - how to interpret results of an experiment or observation
# 4. Numerical methods - solving 1/2/3 when an analytical solution is not possible
# 
# You will use each of these techniques in your future classes! 
# * 06-261 Fluid Mechanics: Navier-Stokes Equations (2)
# * 06-364 Reaction Engineering: e.g. designing reactors (1,2,3)
# * 06-361 Unit Operations: e.g. distillation columns (1,2)
# * 06-363 Heat and Mass Transfer: heat exchangers (1,2)
# * 06-421 Design: solving a process flowsheet (1,2,3)

# # ChemE Motivation for Matrix Algebra
# 
# Consider the following flow diagram:
# 
# <img src="https://docs.google.com/drawings/d/e/2PACX-1vTUqxGdmYDuTmQFyOJVtswSiG9asKjwwk-hBgnalmg430--NVZwZvcq4jGNGWLUgDjz3LctMBG3GMsA/pub?w=1347&amp;h=558">

# So, what do we have here?
# * Three components in each stream
# * Six stream of unknown composition $x_{i=stream}^{j=component}$
# 
# Q: What is the component of each stream?
# 
# Problem: There are 18 compositions and we know only 1: $x_3^c=0$. Therefore, we need 17 equations to solve for 17 unknown compositions!
# 
# Approach: Mole balances! Assume steady state:
# \begin{align*}
# INPUT+GENERATION=OUTPUT
# \end{align*}
# 
# Balances on A:
# \begin{align*}
# Mixer && 95+x_5^A&=x_1^A\\
# Reactor && x_1^A-0.3x_1^A&=x_2^A\\
# Flash && x_2^A&=x_3^A+x_4^A\\
#  && x_3^A&=0.02x_2^A\\
# Split && x_4^A &= x_5^A + x_6^A\\
#  && x_6^A&=0.06x_4^A
# \end{align*}
# 
# Balances on B:
# \begin{align*}
# Mixer && x_5^B&=x_1^B\\
# Reactor && x_1^B+0.3x_1^A&=x_2^B\\
# Flash && x_2^B&=x_3^B+x_4^B\\
#  && x_3^B&=0.98x_2^B\\
# Split && x_4^B&=x_5^B+x_6^B\\
#  && x_6^B&=0.06x_4^B
# \end{align*}
# 
# Balances on C:
# \begin{align*}
# Mixer && 5+x_5^C&=x_1^C\\
# Reactor && x_1^C&=x_2^C \text{  (inert)}\\
# Flash && x_2^C&=x_4^C  \text{  since }x_3^C=0\\
# Split && x_4^C&=x_5^C+x_6^C\\
#  && x_6^C&=0.06x_4^C
# \end{align*}
# 
# * That's a mess! You could solve algebraically by brute force, but we can do this far more quickly and systematically with linear (matrix) algebra! 
# 
# * It's also unclear if we have correctly found all of the necessary equations to solve this system, but we'll be able to show this later!
# 
# * Since all 17 equations are linear (that is, no derivatives since we're at steady state, and only powers of 1)
# 
# * This is quite a simple flow sheet; real chemical plants have about $10^5-10^6$ unknowns and equations!!
# 
# 
# 

# # Matrix/Linear Algebra
# 
# This is an example of a vector:
# \begin{align*}
# \vec{b}=
# \begin{bmatrix}
# 3\\
# 2\\
# -1
# \end{bmatrix}
# \end{align*}
# 
# * The single underline means one dimension
# * The elements are each number in the brackets
# * The length is 3
# * The size is 3x1 (we can think of a vector as a 3x1 array)
# 
# This is an example of a matrix:
# \begin{align*}
# \arr{A}=
# \begin{bmatrix}
# 2&1&3\\
# 1&-2&0
# \end{bmatrix}
# \end{align*}
# 
# * Matrices are denoted by the upper case double underline in notes, and bold in the text
# * Elements are the numbers in the array
# * This is a **2x3 matrix**! Two rows, three columns
# 
# 

# ## Systems of equations as matrices
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
# ## Example
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
# This uses matrix multiplication, which we'll talk about in the next lecture!

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

# In[ ]:




