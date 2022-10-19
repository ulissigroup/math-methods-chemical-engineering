#!/usr/bin/env python
# coding: utf-8

# {{ badge }}
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$

# ## Method of Undetermined Coefficients (sec 2.9)
# \begin{align}
# y'' + ay' + by = r(x)
# \end{align}
# * If $r(x)$ is a function s.t. $r'(x)$ is like $r(x)$, e.g. $e^{\lambda x}$, $\sin(\omega x)$, $\cos(\omega x)$, $x^n$... we can choose $y_P(x)$ based on $r(x)$ \
#  

# ## Rules (from textbook)
# 
# 
# Table 2.1: 
# 
# If $r(x)=$ | Choose $y_P(x)=$
# --- | ---
# $$ke^{rx}$$ | $$ce^{rx}$$
# $kx^n$ ($n=0,1,2...$)| $K_nx^n + K_{n-1}x^{n-1} + ... + K_1x + K_0$
# $k\cos(\omega x)$ or $k\sin(\omega x)$ | $K\cos(\omega x) + M\sin(\omega x)$
# $ke^{\alpha x}\cos(\omega x)$ or $ke^{\alpha x}\sin(\omega x)$ | $e^{\alpha x}(K\cos(\omega x) + M\sin(\omega x))$
# 
# 
# Here are some additional rules; we'll see why these are important later:
# 
# * **Basic Rule**. If $r(x)$ is one of the functions in the first column in Table 2.1, choose in the same line and determine its undetermined coefficients by substituting $y_p$ and its derivatives into the differential equation.
# * **Modification Rule**. If a term in your choice for $y_p$ happens to be a solution of the homogeneous ODE corresponding to (4), multiply this term by x (or by  $x^2$ if this solution corresponds to a double root of the characteristic equation of the homogeneous ODE).
# * **Sum Rule**. If is a sum of functions in the first column of Table 2.1, choose for the sum of the functions in the corresponding lines of the
# second column.

# $\underline{\text{Recap}}$:
# * Given $y'' + ay' + by = r(x) \hspace{3cm} (*)$
# * Solve $y'' + ay' + by = 0$ to get $y_H(x) = c_1y_1 + c_2y_2$
# * Then choose $y_P(x)$ based on $r(x)$
# \begin{align}
# y(x) = c_1y_1 + c_2y_2 + y_P\\
# y'(x) = c_1y_1' + c_2y_2' + y_p'\\
# y''(x) = c_1y_1'' + c_2y_2'' + y_p''
# \end{align}
# * Plug back into $(*)$ to check:
# \begin{align}
# (c_1y_1'' + c_2y_2'' + y_P'') + a(c_1y_1' + c_2y_2' + y_P') + b(c_1y_1 + c_2y_2 + y_P) = r(x)\\
# c_1(y_1'' + ay_1' + by_1) + c_2(y_2'' + ay_2' + by_2) + y_P'' + ay_P' + by_P = r(x)
# \end{align}
# Therefore, $y_P'' + ay_P' + by_P = r(x)$ and $y_P$ is a linear solution to the problem (which is why only some functions work)

# ### Example: $y'' + 2y' - 3y = 4e^{2x}$
# 

# 
# \begin{align}
# y(x) = y_H(x) + y_P(x)
# \end{align}
# 1. Find $y_H(x)$
# \begin{align}
# \lambda^2 + 2\lambda - 3 &= 0\\
# (\lambda + 3)(\lambda - 1) &= 0 \implies \lambda_1 = 1; \hspace{0.5cm} \lambda_2 = -3
# \end{align}
# \begin{align}
# \therefore y_H(x) = c_1e^x + c_2 e^{-3x}
# \end{align}
# 2. Choose $y_P(x)$ based on $r(x)$\
# $\implies$ Since $r(x) = 4e^{2x}$, let's pick $y_p(x) = Ke^{2x}$, where $K$ is the undetermined coefficient (not arbitrary)
# 3. Find $K$ from non-homogeneous ODE
# \begin{align}
# y_P'' + 2y_p' - 3y_p = 4e^{2x}\\
# y_p' = 2Ke^{2x}\\
# y_P'' = 4Ke^{2x}
# \end{align}
# \begin{align}
# \implies 4Ke^{2x} + 4Ke^{2x} - 3Ke^{2x} = 4e^{2x}\\
# 5K = 4\\
# K = \frac{4}{5}
# \end{align}
# 4. Write down general solution
# \begin{align}
# y(x) = y_H(x) + y_P(x)\\
# y(x) = c_1e^x + c_2e^{-3x} + \frac{4}{5}e^{2x}
# \end{align}
# Still need 2 IC's to find $c_1$ and $c_2$

# ## Non-homogeneous Linear $2^\circ$ ODEs Continued

# $\underline{\text{Caution}}$: If you choose a $y_P(x)$ that is not linearly independent of $y_1$ and $y_2$, things will go wrong.
# 
# ### Example: $y'' - 4y = 4e^{2x}$
# 

# 1. Solve $y'' - 4y = 0$
# \begin{align}
# \lambda^2 - 4 &= 0 \implies \lambda_1 = 2; \hspace{0.5cm} \lambda_2 = -2\\
# y_H(x) &= c_1e^{2x} + c_2 e^{-2x}
# \end{align}
# 2. Choose $y_P(x)$ based on $r(x)$\
# Since $r(x) = 4e^{2x}$, choose $y_P(x) = Ke^{2x}$ as before.\
# SPOILER: $Ke^{2x}$ is not linearly independent of $y_1(x) = e^{2x}$
# 3. Find $K$ from non-homogeneous ODE
# \begin{align}
# y_P'' - 4y_P = 4e^{2x}\\
# y_P' = 2Ke^{2x}\\
# y_P'' = 4Ke^{2x}\\
# \implies 4Ke^{2x} - 4Ke^{2x} = 4e^{2x}\\
# 0 = 4e^{2x} \hspace{2cm} ??
# \end{align}
# Step back:
# * If $y_p$ isn't linearly independent, it can be absorbed into $y_H(x)$.
# 
# \begin{align}
# y(x) &= c_1e^{2x} + c_2e^{-2x} + Ke^{2x}\\
# &= (c_1 + K)e^{2x} + c_2 e^{-2x} \implies \text{no undefined coefficient}
# \end{align}
# * Solution: If $y_P(x)$ is linearly dependent on $y_1(x)$ or $y_2(x)$, multiply $y_P(x)$ by $x$ to obtain a linearly independent new $y_P(x)$ based on Reduction of Order.\
# Therefore, choose $y_P(x) = Kxe^{2x}$. Then,
# \begin{align}
# y_P'(x) &= Ke^{2x} + 2Kxe^{2x}\\
# y_P''(x) &= 2Ke^{2x} + 2Ke^{2x} + 4Kxe^{2x}\\
# &= 4Ke^{2x} + 4Kxe^{2x}
# \end{align}
# 
# 3. Find $K$ 
# \begin{align}
# y_P'' - 4y_P = 4e^{2x}\\
# 4Ke^{2x} + 4Kxe^{2x} - 4Kxe^{2x} = 4e^{2x}\\
# 4K = 4\\
# K=1\\
# \therefore y_P(x) = xe^{2x}
# \end{align}
# 4. General Solution
# \begin{align}
# y(x) &= y_H(x) + y_P(x)\\
# y(x) &= c_1e^{2x} + c_2e^{-2x} + xe^{2x}
# \end{align}

# * Another place to be careful:\
# $\implies$ Double roots yields a homogeneous solution
# \begin{align}
# y_H(x) = c_1e^{\lambda x} + c_2xe^{\lambda x}
# \end{align}
# $\implies$ Again, you must be careful to select $y_P(x)$ s.t. it is linearly independent.\
# $\implies$ if $r(x) \propto e^{\lambda x} $ or $xe^{\lambda x}$, then $y_P(x) = Kx^2e^{\lambda x}$ (Obtained from a second reduction of order)

# ### More complex example: $y'' + 2y' -3y = -3x^2 + \sin(x)$
# 
# 

# 
# 1. Solve the homogeneous ODE
# \begin{align}
# y'' + 2y' - 3y = 0\\
# \lambda^2 + 2\lambda -3 = 0\\
# (\lambda+3)(\lambda-1) = 0 && \implies \lambda_1 = 1, \hspace{0.5cm} \lambda_2 = -3\\
# y_H(x) = c_1e^x + c_2e^{-3x}
# \end{align}
# 2. Choose $y_P(x)$ based on $r(x)$\
# $\underline{\text{Rule}}$: If $r(X)$ is a sum of wo functions, then $y_P(x)$ should be a sum of the corresponding functions.
# \begin{align}
# r(x) &= -3x^2 + \sin(x)\\
# \therefore y_P(x) &= Ax^2 + Bx + C + K\sin(x) + M\cos(x)
# \end{align}
# 3. Find undetermined coefficients from ODE
# \begin{align}
# y_P' &= 2Ax + B + K\cos(x) - M\sin(x)\\
# y_P'' &= 2A - K\sin(x) - M\cos(x)
# \end{align}
# * Plugging into ODE:
# \begin{align}
# 2A - K\sin(x) - M\cos(x) + 2[2Ax + B + K\cos(x) - M\sin(x)] - 3[Ax^2 + Bx + C + K\sin(x) + M\cos(x)] = -3x^2 + \sin(x)
# \end{align}
# * Rearrange to collect 'like' terms
# \begin{align}
# -3Ax^2 + (4A - 3B)x + (2A + 2B - 3C) + (-K -2M -3K)\sin(x) + (-M + 2K -3M)\cos(x) = -3x^2 + \sin(x)
# \end{align}
# 
# * Determine the values of the coefficients on each $f(x)$ that will solve the equation
#   1. $-3Ax^2 = -3x^2 \implies A = 1$ 
#   2. $(4A-3B)x - 0$ because $r(x)$ has no '$x$' term . Hence, $4-3B = 0 \implies B = \frac{4}{3}$
#   3. $2A + 2B - 3C = 0 \implies 2+\frac{8}{3} - 3C = 0 \implies C=\frac{14}{9}$
#   4.  $-K -2M -3K = 1 \implies 2M = -4K -1 \implies M = -2K -\frac{1}{2}$
#   5. 
#   \begin{align}
#   -M + 2K -3M = 0\\
#   2K = 4M\\
#   K = 2(-2K-\frac{1}{2}) = -4K-1\\
#   \implies K = -\frac{1}{5}\\
#   \implies M = -\frac{1}{10}
#   \end{align} 
#   \begin{align}
#   \therefore y_P(x) = -3x^2 + \frac{4}{3}x + \frac{14}{9}-\frac{1}{10}(2\sin x + \cos x)
#   \end{align}
# 
# 4. Write down the general solution
# \begin{align}
# y(x) &= y_H(x) + y_P(x)\\
# &= c_1e^x + c_2e^{-3x} + y_P(x)
# \end{align}

# ## Variation of Parameters
# * Non-homogeneous, linear $2^\circ$ ODEs are solvable with Method of Undetermined Coefficients only when $r(x)$ is one of the functions discussed.
# * If $r(x)$ is something different, we need an alternate method $\rightarrow$ variation of parameters.
# * Requires calculation of the Wronskian of the ODE
#   * The value of the Wronskian answers the question: Does a solution exist for the homogeneous equation $y'' + p(x)y' + q(x)y = 0$ ? (specifically at initial condition $x=x_0$)\
# For $y(x) = c_1y_1(x) + c_2y_2(x)$, the Wronskian is defined as a determinant:
# \begin{align}
# W(y_1,y_2) \equiv \left|\begin{array}{} y_1 & y_2 \\ y_1' & y_2'\end{array}\right| = y_1y_2' - y_2y_1'
# \end{align}
# * If $W(y_1,y_2) = 0$ at $x=x_0$, then $y_1$ and $y_2$ are linearly dependent $\rightarrow$ not good
# * If $W(y_1,y_2) \neq 0$ on the interval of interest, which includes $x_0$, then the solutions $y_1$ and $y_2$ are linearly independent.
# 
# 
# 

# ### Example
# 

# Constant coefficient case with 2 distinct roots: $y_1=\lambda_1e^{\lambda_1 x}$ and $y_2 = \lambda_2e^{\lambda_2 x}$
# \begin{align}
# \therefore y_1' = \lambda_1 e^{\lambda_1x} && \& && y_2' = \lambda_2e^{\lambda_2x}
# \end{align}
# \begin{align}
# W(y_1,y_2) &= e^{\lambda_1 x}\cdot \lambda_2e^{\lambda_2x} - e^{\lambda_2 x} \cdot \lambda_1 e^{\lambda_1x}\\
# &= (\lambda_2 - \lambda_1)e^{(\lambda_1 + \lambda_2)x} 
# \end{align}
# which never equals zero except for $\lambda_2 = \lambda_1$ which by definition is not the case. Therefore, solution exists for all initial conditions.

# ## Variation of Parameters (sec 2.10)
# For solving $y'' + p(x)y' + q(x)y = r(x)$\
# where $p(x)$, $q(x)$ and $r(x)$ are continuous over some interval $I$.\
# $\implies$ Based on same idea as reduction of order
# * First, find $y_H(x) = c_1y_1 + c_2y_2$
# * Then we will base $y_P(x)$ on $y_H(x)$, $\underline{\text{not}}$ $r(x)$ like before
# \begin{align}
# y_P(x) = u(x)y_1(x) + v(x)y_2(x) \hspace{3cm} (*)
# \end{align}
# where $u(x)$ and $v(x)$ are "parameters" that replace $c_1$ and $c_2$ in $y_H(x)$
# \begin{align}
# y_P' = u'y_1 + uy_1' + v'y_2 + vy_2'
# \end{align}
#   * Since there is only one equation $(*)$ for $y_P(x)$ and two unknown functions ($u(x)$ and $v(x)$), there will likely be many choices of $u$ and $v$ that will work.
#   * However, we can impose an additional constraint (since there is an extra degree of freedom) if we like
#   * Choose this constraint to help simplify the solution method
#   \begin{align}
#   u'y_1 + v'y_2 = 0
#   \end{align}
#   That simplifies our first derivative
#   \begin{align}
#   y_p' = uy_1' + vy_2'
#   \end{align}
#   then,
#   \begin{align}
#   y_p'' = uy_1'' + u'y_1' + vy_2'' + v'y_2'
#   \end{align}
# * We can then plug these terms back into non-homogeneous ODE:
# \begin{align}
# y_p'' + p(x)y_p' + q(x)y_p = r(x) \\
# (uy_1'' + u'y_1' + vy_2'' + v'y_2') + p(x)(uy_1' + vy_2') + q(x)(uy_1 + vy_2) = r(x)
# \end{align}
# * Rearrange strategically:
# \begin{align}
# u(y_1'' + py_1' + qy_1) + v(y_2'' + py_2' + qy_2) + u'y_1' + v'y_2' = r(x)
# \end{align}
# By definition, $y_1$ and $y_2$ are solutions of homogeneous equation so $y_1''+py_1'+qy_1=0$ and $y_2''+py_2'+qy_2$. Thus
# \begin{align}
# \therefore u'y_1' + v'y_2' &= r\\
# u'y_1 + v'y_2 &= 0 \text{ [from above]}
# \end{align}
# We now have two equations with two unknowns $u$ and $v$
# \begin{align}
# \begin{bmatrix}y_1 & y_2\\ y_1' & y_2'\end{bmatrix} \begin{bmatrix} u' \\ v' \end{bmatrix} = \begin{bmatrix} 0 \\ r \end{bmatrix}
# \end{align}
# * Could solve by Gauss Elimination or by something called Cramer's rule, which says\
# for $\begin{bmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} b_1 \\b_2 \end{bmatrix}$\
# that
# \begin{align}
# x_1 =\frac{ \left| \begin{array}{} b_1 & a_{12} \\ b_2 & a_{22}\end{array}\right|} 
# {\left| \begin{array}{} a_{11} & a_{12} \\ a_{21} & a_{22}\end{array}\right| } \hspace{0.5cm}\text{&} \hspace{0.5cm}
# x_2 = \frac{ \left| \begin{array}{} a_{11} & b_1 \\ a_{21} & b_2\end{array}\right|} 
# {\left| \begin{array}{} a_{11} & a_{12} \\ a_{21} & a_{22}\end{array}\right| }
# \end{align}
# 
# 
# 

# \begin{align}
# \therefore u' = \frac{ \left| \begin{array}{} 0 & y_2 \\ r & y_2'\end{array}\right|} 
# {\left| \begin{array}{} y_1 & y_2 \\ y_1' & y_2'\end{array}\right| }
# = \frac{-ry_2}{W}
# \end{align}
# and 
# \begin{align}
# v' = \frac{ \left| \begin{array}{} y_1 & 0 \\ y_1' & r\end{array}\right|} 
# {\left| \begin{array}{} y_1 & y_2 \\ y_1' & y_2'\end{array}\right| }
# = \frac{ry_1}{W}\\
# \therefore u(x) = \int \frac{-r(x)y_2(x)}{W}dx\\
# v(x) = \int \frac{r(x)y_1(x)}{W}dx
# \end{align}
# and then
# \begin{align}
# y(x) &= y_H(x) + y_P(x)\\
# y(x) &= c_1y_1 + c_2y_2 + uy_1 + vy_2
# \end{align}
# Therefore, V.O.P can be written as:
# \begin{align}
# y(x) = \left[c_1 + \int \frac{-ry_2}{W} dx \right]y_1 + \left[c_2 + \int \frac{ry_1}{W} dx \right] y_2 && \text{general solution}
# \end{align}
# * Nice solution technique because it's general and can be used for any $r(x)$

# ### Example 
# 
# First, let's do one that can be solved with MUC or VoP
# \begin{align}
# y'' + y = x
# \end{align}
# 
# 

# for either solution method, solve $y''+y=0$ first
# \begin{align}
# \lambda^2 + 1 = 0 \implies \lambda^2 = -1 \implies \lambda = \pm i\\
# a = 0, \hspace{0.5cm} b = 1 \implies \omega = \sqrt{b - \frac{1}{4}a^2} = 1\\
# \end{align}
# \begin{align}
# y_H(x) &= e^{-\frac{a}{2}}[c_1\sin(\omega x) + c_2\cos(\omega x)]\\
# &= c_1\sin x + c_2 \cos x
# \end{align}
# 1. if using MUC, choose
# \begin{align}
# y_P(x) = Ax + B\\
# y_P' = A\\
# y_P'' = 0\\
# \end{align}
# substitute,
# \begin{align}
# y_P'' + y_P = x\\
# 0 + Ax + B = x \\
# \implies A = 1, \hspace{0.5cm} B = 0\\
# \therefore y_P(x) = +x
# \end{align}
# and 
# \begin{align}
# y(x) = c_1\sin x + c_2\cos x + x
# \end{align}
# 2. if using VoP, then 
# \begin{align}
# y_P(x) = u(x)\sin(x) + v(x)\cos(x)
# \end{align}
# where
# \begin{align}
# u(x) = \int \frac{-r(x)y_2(x)}{W}dx\\
# v(x) = \int \frac{r(x)y_1(x)}{W}dx
# \end{align}
# * We know $r(x)$, $y(x)$ and $y_2(x)$. Let's find the Wronskian
# \begin{align}
# y_1 = \sin x\\
# y_1' = \cos x \\
# y_2 = \cos x\\
# y_2' = -\sin x
# \end{align}
# Hence, 
# \begin{align}
# W(y_1, y_2) &= y_1y_2' - y_2y_1'\\
# &= -\sin^2x - \cos^2x\\
# &= -(\sin^2x + \cos^2x)\\
# &= -1\\
# \therefore u(x) &= \int r(x)y_2(x)dx\\
# &= \int x\cos x dx\\
# &=\cos x + x\sin x\\
# v(x) &= -\int r(x) y_1(x) dx\\
# &=-\int x \cdot \sin x dx\\
# &= -\sin x + x\cos x && \text{(integration by parts)}
# \end{align}
# 
# Then, 
# \begin{align}
# y_P &= u\sin x + v\cos x\\
# &= (\cos x + x\sin x)\cdot \sin x + (x \cos x - \sin x) \cos x\\
# &= x (\sin^2x + \cos^2x) \\
# &= x
# \end{align}
# and 
# \begin{align}
# \therefore y(x) &= y_H(x) + y_P(x)\\
# y(x) &= c_1\sin(x) + c_2\cos(x) + x && \text{same answer :)}
# \end{align}
# $\underline{\text{Ques}}$: Which method would you prefer in this case? What if we had many undefined coefficients?

# ### Example
# 
# Now an example for which we can't use MUC
# \begin{align}
# y'' -2y' + y = \frac{12e^x}{x^3}
# \end{align}
# 

# 1. Solve homogeneous case
# \begin{align}
# \lambda^2 - 2\lambda + 1 = 0\\
# (\lambda - 1)^2 = 0 \implies \lambda = 1\\
# \therefore y_H(x) = c_1e^x + c_2xe^x
# \end{align}
# 2. Assume $y_P(x) = u(x)y_1(x) + v(x)y_2(x)$
#   * Find Wronskian
#   \begin{align}
#   y_1 = e^x, &\hspace{0.5cm} y_1' = e^x\\
#   y_2 = xe^x, &\hspace{0.5cm} y_2' = e^x(x+1)  
#   \end{align}
#   Hence,
#   \begin{align}
#   W(y_1,y_2) &= y_1y_2' - y_2y_1'\\
#   &= e^{2x}(x+1) - xe^{2x}\\
#   &= xe^{2x} + e^{2x} - xe^{2x}\\
#   W &= e^{2x}
#   \end{align}
#   * Then
#   \begin{align}
#   u(x) &= \int \frac{-r(x)y_2(x)}{W} dx\\
#   &= -\int \frac{12e^x}{x^3} \cdot \frac{xe^x}{e^{2x}}dx\\
#   &= -12 \int x^{-2}dx\\
#   u(x) &= 12x^{-1}\\
#   v(x) &= \int\frac{r(x)y_1(x)}{W} dx\\
#   &= \int \frac{12e^x}{x^3} \cdot \frac{e^x}{e^{2x}} dx\\
#   &= 12\int x^{-3}\\
#   v(x) &= -6x^{-2}
#   \end{align}
# 3. Write down general solution
# \begin{align}
# y(x) &= y_H + y_P\\
# &= (c_1 + u)y_1 + (c_2 + v)y_2\\
# &= \left(c_1 + \frac{12}{x}\right)e^x + \left(c_2 - \frac{6}{x^2}\right)xe^x\\
# y(x) &= (c_1 + c_2x)e^x + \frac{6}{x}e^x
# \end{align}
# Still need two IC's to find particular solution
