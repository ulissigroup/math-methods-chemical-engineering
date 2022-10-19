#!/usr/bin/env python
# coding: utf-8

# {{ badge }}
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$

# # Recap 
# 
# We did a lot last class!
# * Definition of differential equations and classifications
#   * Linear/non-linear
#   * ODE/PDE
#   * Order
#   * Separable diff eq
# * Implicit vs explicit solutions
# * General vs particular solutions
# * Mol balance applications with an accumulation term
# 
# Misc logistics:
# * Exam practice is up on google drive (Exam Thursday next week, 2/24)
# * First SI session will be this week (I think weekend?) and cover a review for exam 1
# * There will be a shorter HW assigned this week to practice exact differential equations

# # Exact $1^\circ$ ODEs
# $\rightarrow$ Specific type of ODE\
# $\rightarrow$ Method of solving exact $1^\circ$ ODEs was developed by guessing and messing around.\
# $\rightarrow$ Now, if we can classify as "exact", we have a method to find the solution.
# * The $1^\circ$ ODE $\frac{dy}{dx} = f(x,y)$ is exact if:
#   1. It can be written as $M(x,y)dx + N(x,y)dy = 0 \hspace{3cm} (*)$
#   2. There is some function $u(x,y)$ such that 
#   \begin{align}
#   \frac{\partial u}{\partial x} =M(x,y) && \text{and} && \frac{\partial u}{\partial y}=N(x,y)
#   \end{align}
# $\rightarrow$ plugging this into $(*)$, we see:
# \begin{align}
# \frac{\partial u}{\partial x}dx + \frac{\partial u}{\partial y} dy = 0\\
# \therefore du = 0
# \end{align}
# $\rightarrow$ integrating this gives $u(x,y) = c$\
# $\rightarrow$ Therefore if we have an exact equation $y'=f(x,y)$, then the general solution can be written as $u(x,y) = c \implies$ function of $x,y$ gives a constant, not a derivative.
# 
# To prove exactness:
# * Does $u(x,y)$ exist such that 
# \begin{align}
# \frac{\partial u}{\partial x} = M(x,y) && \text{and} && \frac{\partial u}{\partial y} = N(x,y) \hspace{1cm} ?
# \end{align}
# 
# This strategy can also be applied to higher order differential equations, but it's not as helpful there so I'm not going to cover that in lectures.

# * If we assume continuity of the derivatives of $u(x,y)$, then we can write
# \begin{align}
# \frac{\partial^2u}{\partial x\partial y} &= \frac{\partial^2u}{\partial y\partial x}\\
# \frac{\partial}{\partial x}\left(\frac{\partial u}{\partial y} \right) &= \frac{\partial}{\partial y}\left(\frac{\partial u}{\partial x} \right)\\
# \frac{\partial}{\partial x}N(x,y) &= \frac{\partial}{\partial y}M(x,y) && \leftarrow \text{if exact, this will be true} 
# \end{align}
# 

# ### Example
# $$3xy^2y' + y^3 = 0$$

# 
# $3xy^2y' + y^3 = 0 \rightarrow$ we can solve 2 ways: by separating and by exactness
# 1. Separate
# \begin{align}
# 3xy^2 \frac{dy}{dx} = -y^3\\
# \frac{3}{y} dy = -\frac{1}{x}dx\\
# 3 \ln y = - \ln x + \beta\\
# \ln (y^3) + \ln x = \beta\\
# \ln(xy^3) = \beta\\
# xy^3 = \exp \beta \\
# xy^3 = \alpha &&\rightarrow \text{implicit solution(didn't solve for y) by separation}
# \end{align}
# 
# 2. Let's check to see if we can solve by exactness:\
#   a) Write as $M(x,y)dx + N(x,y)dy=0$:
#   \begin{align}
#   y^3dx + 3xy^2 dy = 0
#   \end{align}
#   b) Check if $\frac{\partial}{\partial x}N(x,y) = \frac{\partial}{\partial y} M(x,y)$
#   \begin{align}
#   \frac{\partial}{\partial x} (3xy^2) &= \frac{\partial}{\partial y}(y^3)\\
#   3y^2 &= 3y^2 && \rightarrow \text{ODE is exact; we can continue}
#   \end{align}
#   c) Find $u(x,y)$
#   \begin{align}
#   \text{we know that} && \frac{\partial u}{\partial x} = M(x,y) = y^3\\
#   \text{and} && \frac{\partial u}{\partial y} = N(x,y) = 3xy^2\\
#   \end{align}
#   Since these are partial derivatives, we can integrate them as follows:
#   \begin{align}
#   \frac{\partial u}{\partial x} &= y^3\\
#   \int du &= \int y^3 dx\\
#   u(x,y) &= xy^3 + k(y)
#   \end{align}
#   we need to include the last term because its derivative wr.t. $x=0$\
#   We can find $k(y)$ by considering our second partial derivative:
#   \begin{align}
#   \frac{\partial u}{\partial y} = 3xy^2\\
#   \frac{\partial}{\partial y}(xy^3 + k(y)) = 3xy^2\\
#   3xy^2 + \frac{dk}{dy} = 3xy^2\\
#   \frac{dk}{dy} = 0 \implies k(y) = k && \text{(constant)}\\
#   \therefore u(x,y) = xy^3 + k
#   \end{align}
#   d) Use $u(x,y)$ to identify solution
#   \begin{align}
#   u(x,y) = c = xy^3 + k \\
#   \implies xy^3 = \alpha
#   \end{align}
#   Implicit solution; same as when we separated.

# ### Example 
# 
# \begin{align}
# y' = \frac{-2xy^3 - 2}{3x^2y^2 + e^y} && \text{(not separable)}
# \end{align}

# 
# 1. Put into the form $M(x,y)dx + N(x,y)dy = 0$
# \begin{align}
# (2xy^3 + 2) dx + (3x^2y^2 + e^y)dy = 0
# \end{align}
# 2. Check for exactness
# \begin{align}
# \frac{\partial}{\partial x}N(x,y) &= \frac{\partial}{\partial y}M(x,y) && \text{(must be true)}\\
# \frac{\partial}{\partial x} (3x^2y^2 + e^y) &= \frac{\partial}{\partial y} (2xy^3 + 2)\\
# 3y^2\frac{\partial}{\partial x} (x^2) &= 2x \frac{\partial}{\partial y}(y^3)\\
# 6y^2x &=6xy^2 &&\text{(exact)}
# \end{align}
# 3. Now find $u(x,y)$ such that 
# \begin{align}
# \frac{\partial u}{\partial y} = N(x,y)=3x^2y^2 + e^y && \text{and} && \frac{\partial u}{\partial x} = M(x,y)=2xy^3 + 2 && \text{(can start with either)}
# \end{align}

# \begin{align}
# \int \partial u = \int (3x^2y^2 + e^y) dy\\
# u(x,y) = x^2y^3 + e^y + k(x)
# \end{align}
# $\rightarrow$ determine $k(y)$ with second partial derivative:
# \begin{align}
# \frac{\partial}{\partial x}(x^2y^3 + e^y + k(x)) &= 2xy^3 + 2\\
# 2xy^3 + \frac{dk}{dx} &= 2xy^3 + 2\\
# \frac{dk}{dx} &= 2\\
# k(x) &= 2x + k
# \end{align}
# $\therefore u(x,y) = x^2y^3 + e^y + 2x + k$
# 4. Find implicit solution using $u(x,y) = c$
# \begin{align}
# c = x^2y^3 + e^y + 2x + k\\
# x^2 y^3 + e^y + 2x = \alpha && \text {general solution}
# \end{align}
# 

# # Integrating Factors
# * Intefrating Factors (IFs) are a tool we can use to make some inexact equations exact so we can solve.\
# In general, we have
# \begin{align}
# P(x,y)dx + Q(x,y) dy = 0
# \end{align}
# Where $\frac{\partial P}{\partial y} \neq \frac{\partial Q}{\partial x} \implies$ not exact
# * We can try to find an I.F., $F(x,y)$, s.t. $F(x,y)P(x,y)dx + F(x,y)Q(x,y)dy = 0$ is an exact equation. This would require:
# \begin{align}
# \frac{\partial}{\partial y}[F(x,y)P(x,y)] &= \frac{\partial}{\partial x} [F(x,y)Q(x,y)]\\
# F(x,y)\frac{\partial P}{\partial y} + P(x,y)\frac{\partial F}{\partial y} &= F(x,y)\frac{\partial Q}{\partial x} + Q(x,y)\frac{\partial F}{\partial x}
# \end{align}
# Knowing only P(x,y) and Q(x,y) would make this nearly impossible to solve. Can use the following strategy:
# 1. First assume that $F(x,y) = F(x)$ only. Now, 
# \begin{align}
# \frac{\partial F}{\partial y} = 0 && \text{and} && \frac{\partial F}{\partial x} = \frac{dF}{dx}  && \text{(partial becomes full)}
# \end{align}
# \begin{align}
# \therefore F(x) \frac{\partial P}{\partial y} + P(x,y) \frac{\partial F}{\partial y} &= F(x) \frac{\partial Q}{\partial x} + Q(x,y) \frac{dF}{dx}\\
# F(x)\left[\frac{\partial P}{\partial y} - \frac{\partial Q}{\partial x} \right] &= Q(x,y)\frac{dF}{dx}\\
# \frac{1}{F(x)}\frac{dF}{dx} &= \frac{1}{Q(x,y)}\left[\frac{\partial P}{\partial y} - \frac{\partial Q}{\partial x} \right]
# \end{align}
# we know the RHS term, call it $R(x,y)$\
# $\rightarrow$ If $R(x,y) = R(x)$ only, then you can find $F(x)$ from integration
# \begin{align}
# \frac{\partial F}{F} = R(x)dx \implies F(x) = \exp\left[\int R(x)dx \right]
# \end{align}
# $\rightarrow$ If $R(x,y)\neq R(x)$, then out assumption that $F(x,y)=F(x)$ is not true. Then try step 2.

# 2. Assume $F(x,y) = F(y)$ only. Then,
# \begin{align}
# \frac{\partial F}{\partial x} = 0 && \text{and} && \frac{\partial F}{\partial y} = \frac{dF}{dy}  
# \end{align}
# Partial derivative equation simplifies to:
# \begin{align}
# \frac{1}{F} \frac{dF}{dy} = \frac{1}{P(x,y)}\left[\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right]
# \end{align}
# Again, let RHS = $R(x,y)$\
# $\rightarrow$ If $R(x,y) = R(y)$ only then $F(y) = \exp \left[\int R(y)dy \right]$\
# $\rightarrow$If $R(x,y) \neq R(y)$, then $F(x,y)$ = function of both $x$ and $y$ and you will need some other clue to solve.

# ### Example
# 
# \begin{align}
# \frac{dy}{dx} = \frac{-2xy}{4y + 3x^2}
# \end{align}

# 
# 1. rewrite as: $(2xy)dx + (4y + 3x^2) dy = 0$
# 2. check for exactness:
# \begin{align}
# \frac{\partial P}{\partial y} &= \frac{\partial Q}{\partial x}\\
# 2x &\neq 6x \implies \text{not exact}
# \end{align}
# 3. Try to find an I.F. $F(x,y)=F(x)$
# \begin{align}
# \frac{\partial (FP)}{\partial y} &= \frac{\partial (FQ)}{\partial x}\\
# F\frac{\partial P}{\partial y} + P\frac{\partial F}{\partial y} &= F \frac{\partial Q}{\partial x} + Q\frac{\partial F}{\partial x}\\
# \therefore \frac{1}{F} \frac{dF}{dx} &= \frac{1}{Q}\left(\frac{\partial P}{\partial y} - \ \frac{\partial Q}{\partial x}\right)\\
# &=\frac{1}{4y + 3x^2}(2x-6x)\\
# &=-\frac{4x}{4y+3x^2} \implies \text{($R(x,y)$ depends on $x$ and $y\rightarrow$ no good)}
# \end{align}
# 4. Try to find I.F. $F(x,y) = F(y)$
# \begin{align}
# \frac{1}{F} \frac{dF}{dy} &= \frac{1}{P}\left(\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} \right)\\
# &=\frac{1}{2xy}(6x-2x) = \frac{2}{y} = R(y) && \text{(good!)}
# \end{align}
# Continuing,
# \begin{align}
# \int \frac{dF}{F} &= \int \frac{2}{y}dy\\
# \ln F(y) &= 2 \ln y\\
# F(y) &= y^2
# \end{align}
# 5. We can now make our equation exact using $F(y)$
# \begin{align}
# F(y)[P(x,y)dx + Q(x,y)dy] = 0\\
# y^2 [(2xy)dx + (4y + 3x^2)dy] = 0\\
# 2xy^3 dx + (4y^3 + 3x^2y^2)dy = 0
# \end{align}
# 6. Check for exactness:
# \begin{align}
# \frac{\partial}{\partial y}(2xy^3) &= \frac{\partial}{\partial x}(4y^3 + 3x^2y^2)\\
# 6xy^2 &=6xy^2 && \text{exact}
# \end{align}
# 7. Find $u(x,y)$
# \begin{align}
# \frac{\partial u}{\partial x} = 2xy^3 \implies u(x,y) = x^2y^3 + k(y)\\
# \frac{\partial u}{\partial y} = \frac{\partial}{\partial y} (x^2y^3 + k(y)) = 4y^3 + 3x^2y^2\\
# 3x^2y^2 + \frac{dk}{dy} = 4y^3 + 3x^2y^2\\
# k(y) = y^4 + k
# \end{align}
# 8. Find implicit solution using $u(x,y) = c$ 
# \begin{align}
# x^2y^3 + y^4 + k = c\\
# x^2 y^3 + y^4 = \alpha && \text{general solution}
# \end{align}

# In[ ]:




