#!/usr/bin/env python
# coding: utf-8

# {{ badge }}
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$

# # Recap from last class 
# 
# For the 2nd order, linear, constant coefficients, homogeneous differential equation $y'' + ay' + by = 0$ where $a^2-4b<0$, the general solution is:
# \begin{align*}
# y(x) &= e^{-\frac{a}{2}x}[A\sin(\omega x) + B\cos(\omega x)] && \omega = \sqrt{b - \frac{1}{4}a^2}
# \end{align*}
# where $A$ and $B$ are our arbitrary constants. All components of the solution are real.
# 

# ## Example
# 
# $y'' + 2y' + 3y = 0$ ; $\hspace{0.5cm}y(0) = 2$ and $y'(0) = -3$\
# 
# 

# Characteristic equation:
# \begin{align}
# \lambda^2 + 2\lambda + 3 = 0\\
# a = 2, \hspace{1cm} b = 3
# \end{align}
# * Since $a^2 - 4b = -8$, we have complex roots
# \begin{align}
# \lambda &= -\frac{1}{2}a \pm i\sqrt{b - \frac{1}{4}a^2}\\
# &= -\frac{1}{2}a \pm i \omega \implies \omega = \sqrt{2}\\
# \lambda &= -1 \pm i\sqrt{2}
# \end{align}
# * Great, so we know our general solution is:
# \begin{align}
# y(x) &= e^{-\frac{a}{2}x} [A\cos(\omega x) + B\sin(\omega x)]\\
# &= e^{-x}[A\cos(\sqrt{2}x) + B \sin(\sqrt{2}x)]
# \end{align}
# * To find $A$ and $B$, we will also need an expression for $y'(x)$
# \begin{align}
# y'(x) &= -e^{-x}[A\cos(\sqrt{2}x) + B\sin(\sqrt{2}x)] + + e^{-x}[-\sqrt{2} A\sin(\sqrt{2}x) + \sqrt{2} B\cos(\sqrt{2}x)]\\
# &= e^{-x}[(\sqrt{2} B - A) \cos(\sqrt{2}x) - (\sqrt{2}A + B)\sin(\sqrt{2}x)]
# \end{align}
# 
# Using Initial Conditions:
# \begin{align}
# y(0) = 2 &= e^{-0}[A\cos(0) + B\sin(0)]\\
# 2&=A
# \end{align}
# and,
# \begin{align}
# y'(0) = -3 &= e^{-0}[(\sqrt{2}B - A) \cos(0) - (\sqrt{2}A + B) \sin(0)]\\
# -3&=\sqrt{2}B - A\\
# B &=\frac{1}{\sqrt{2}}(-3 + A) = \frac{1}{\sqrt{2}}(-3+2)\\
# B &= -\frac{1}{\sqrt{2}} = -\frac{\sqrt{2}}{2}
# \end{align}
# \begin{align}
# \therefore y(x) = e^{-x}[2 \cos(\sqrt{2}x) -\frac{\sqrt{2}}{2}\sin(\sqrt{2}x)] \implies \text{particular solution}
# \end{align}

# $\underline{\text{Summary}}$: Homogeneous, linear $2^\circ$ ODEs with constant coefficients
# \begin{align}
# y'' + ay' + by = 0
# \end{align}
# * Three types of solutions, arising from three cases for the characteristic equation $\lambda^2 + a\lambda + b = 0$
#   1. Two real and different roots
#   \begin{align}
#   y(x) = c_1 ^{\lambda_1x} + c_2 e^{\lambda_2x} && \text{where } \lambda = \frac{1}{2}(-a \pm \sqrt{a^2 - 4b})
#   \end{align}
#   2. Real, double roots
#   \begin{align}
#   y(x) = c_1 e^{\lambda x} + c_2 x e^{\lambda x}
#   \end{align}
#   3. Complex roots
#   \begin{align}
#   y(x) = c_1 e^{-\frac{a}{2}x}[A \sin(\omega x) + B \cos(\omega x)] && \text{where } \omega = \sqrt{b - \frac{1}{4}a^2}
#   \end{align}
# 
# In each case, we have a basis of $\underline{\text{two}}$ linearly independent solutions (for $2^{nd}$ order). 
# $\rightarrow$ in case there's just one, this necessitates that $\lambda_1 \neq \lambda_2$, which is the definition of the case

# ## Non-homogeneous Linear $2^\circ$ ODEs (sec 2.8)
# \begin{align}
# y'' + p(x) y' + q(x) y = r(x)
# \end{align}
# * We build upon our knowledge of solving the homogenous case to solve non-homogeneous case.\
# $\underline{\text{General solution}}$:
# \begin{align}
# y(x) = y_H(x) + y_P(x)
# \end{align}
# where $y_H(x) = c_1y_1(x) + c_2y_2(x)$, which is the "homogeneous solution" for $r(x)\equiv 0$ and $y_P(x)$ is the "particular solution" which accounts for non-homogeneous part and has no arbitrary constants (we need only two for a $2^\circ$ ODE)
# * $\underline{\text{NOTE}}$: $c_1$ and $c_2$ must be found using $y(x) = y_H(x) + y_P(x)$, not just $y_H(x)$
# * We will do an example first and then describe the general method to find $y_P(x)$

# ### Example
#  $y'' + 4y = e^x $ ; $\hspace{0.5cm}y(0) = \frac{3}{5} $ and $y'(0)=\frac{4}{5}$
# 
# 

# 
# * First, solve the homogeneous case:
# \begin{align}
# y'' + 4y = 0 \implies \text{assume solution of form $e^{\lambda x}$}
# \end{align}
# characteristic equation:
# \begin{align}
# \lambda ^ 2 + 4 = 0 \implies \lambda = \pm 2i\\
# a = 0, b = 4 \implies \omega = \sqrt{b - \frac{1}{4}a^2} = 2
# \end{align}
# $\implies$ for imaginary roots, we know the solution is:
# \begin{align}
# y_H(x) &= [c_1 \sin(\omega x) + c_2 \cos(\omega x)] e^{-\frac{a}{2}x}\\
# &= c_1 \sin(2x) + c_2 \cos(2x) \implies \text{solution to homogeneous equation} 
# \end{align}
# Now, I will tell you that
# \begin{align}
# y_P(x) &= \frac{1}{5}e^x\\
# \therefore y(x) &= y_H(x) + \frac{1}{5}e^x\\
# y'(x) &= y_H'(x) + \frac{1}{5}e^x\\
# y''(x) &= y_H''(x) + \frac{1}{5}e^x
# \end{align}
# Plugging these into original non-homogeneous equation:
# \begin{align}
# y_H'' + \frac{1}{5}e^x + 4(y_H + \frac{1}{5}e^x) = e^x\\
# y_H'' + \frac{1}{5}e^x + 4y_H + \frac{4}{5}e^x = e^x\\
# y_H'' + 4y_H + e^x = e^x && \implies \text{equation holds and } y_P(x) = \frac{1}{5}e^x \text{is correct}
# \end{align}
# 
# * Now we need the particular solution to the general solution of non-homogeneous equation using initial conditions.
#   1. 
#   \begin{align}
#   y(0) = \frac{3}{5} &= c_1 \sin(2\cdot 0) + c_2 \cos (2\cdot 0) + \frac{1}{5} e^0\\
#   \frac{3}{5} &= c_2 + \frac{1}{5}\\
#   c_2 &= \frac{2}{5}
#   \end{align}
#   2. We need expression for $y'(x)$ for $2^{nd}$ I.C.
#   \begin{align}
#   y'(x) &= 2c_1\cos(2x) - 2c_2 \sin(2x) + \frac{1}{5}e^x\\
#   y'(0) = \frac{4}{5} &= 2c_1 \cos(2\cdot 0) - 2c_2 \sin(2\cdot 0) + \frac{1}{5}e^0\\
#   \frac{4}{5} &= 2c_1 + \frac{1}{5}\\
#   c_1 &= \frac{3}{10}
#   \end{align}
#   \begin{align}
#   \therefore y(x) = \frac{3}{10} \sin(2x) + \frac{2}{5} \cos(2x) + \frac{1}{5} e^x \implies \text{particular solution to non-homogeneous linear $2^\circ$ ODE}
#   \end{align}
# 
# * How do we determine $y_P(x)$?
# 

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

# # In-class problem
# 
# Consider the differential equation $$y''+y=0.001x^2$$ The homogeneous solution is $$y_h=A\cos x+B\sin x$$, as we saw in class last week. Find a particular solution to this differential equation using the MoUC rules above.

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
# * If $y(p)$ isn't linearly independent, it can't be absorbed into $y_H(x)$.
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
