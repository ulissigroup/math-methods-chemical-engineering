#!/usr/bin/env python
# coding: utf-8

# {{ badge }}
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$

# # Integrating Factors - Recap
# $\rightarrow$ Functions that make inexact $1^\circ$ ODEs exact (can't always find them).\
# $\rightarrow$ Two possibilities in this class: $F(x)$ or $F(y)$
# * if $F(x)$, then
# \begin{align}
# \frac{1}{F}\frac{dF}{dx} = \frac{1}{Q}\left[\frac{\partial P}{\partial y} - \frac{\partial Q}{\partial x}\right] = R(x)
# \end{align}
# and $F(x) = \exp[\int R(x)dx]$
# * if $F(y)$, then 
# \begin{align}
# \frac{1}{F}\frac{dF}{dy} = \frac{1}{P}\left[\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y}\right] = R(y)
# \end{align}
# and $F(y) = \exp[\int R(y)dy]$
# 
# $\rightarrow$ Once you find your I.F., follow the steps for solving an exact equation:
# 1. Write in form $M(x)dx + N(y)dy = 0$
# 2. Check for exactness
# 3. If exact find $u(x,y)$
# 4. Find implicit general solution from $u(x,y)=c$

# # Linear $1^\circ$ ODEs
# A $1^\circ$ differential equation is linear if it can be written:
# \begin{align}
# y' + p(x) y = r(x)
# \end{align}
# $\rightarrow$ it is linear in $y$ and $y'$\
# $\rightarrow$ $p$ and $r$ may be any function of $x$ only
# * Two cases to solve:
# 1. Homogeneous Linear ODE
# \begin{align}
# r(x) = 0 && \implies \text{can solve by separating}
# \end{align}
# 2. Non-homogeneous Linear ODE
# \begin{align}
# r(x) \neq 0 && \implies \text{can solve as an exact equation with a simple integrating factor}
# \end{align}

# ## Homogeneous Case
# \begin{align}
# y' &+ p(x)y = 0\\
# \frac{dy}{dx} &= -p(x) y\\
# \frac{dy}{y} &= -p(x)dx\\
# \ln y &= -\int p(x)dx + \beta\\
# y(x) &= \alpha \exp\left[- \int p(x) dx\right]
# \end{align}
# 

# ## Non-Homogeneous Case
# \begin{align}
# y' &+ p(x)y = r(x)\\
# \frac{dy}{dx} &= -p(x)y + r(x)\\
# \end{align}
# 1. Write in the form
# \begin{align}
# [p(x)y - r(x) ]dx + 1 dy = 0
# \end{align}
# 2. Check for exactness:
# \begin{align}
# \frac{\partial}{\partial y}P(x,y) &= \frac{\partial}{\partial x} Q(x,y)\\
# p(x) & \neq 0 && \implies \text{not exact}
# \end{align}
# 3. Find an integrating factor. Try $F(x,y)=F(x)$
# \begin{align}
# \frac{1}{F}\frac{dF(x)}{dx} &= \frac{1}{Q} \left[\frac{\partial P}{\partial y} - \frac{\partial Q}{\partial x} \right]\\
# &= \frac{1}{1}(p(x) - 0)\\
# &= p(x)\\
# \therefore F(x) &= \exp\left(\int p(x)dx\right) && \implies \text{this will always be the I.F. for a linear ODE}
# \end{align}
# 
# 

# * Once you recognize that you have a linear $1^\circ$ ODE, you know that you can solve it as an exact equation using the I.F. $F(x) = \exp\left( \int p(x)dx \right)$
# * Therefore since we will always have the same I.F., let's solve the non-homogeneous linear $1^\circ$ ODE to obtain a general solution $\implies$ won't have to derive every time.
# \begin{align}
# F(x)y' + F(x) p(x) y = F(x) r(x)\\
# \exp\left[ \int p(x) dx\right] y' + \exp \left[ \int p(x) dx\right] p(x) y = \exp\left[\int p(x)dx \right] r(x)\\
# \frac{d}{dx} \left( \exp\left[\int p(x) dx \right] y \right) = \exp\left[\int p(x)dx \right] r(x) && \text{(chain rule)}
# \end{align}
# * Note that the RHS is dependent only on x, so we can simply integrate 
# \begin{align}
# \exp \left[\int p(x)dx \right] y = \int \exp\left[\int p(x)dx \right] r(x) dx + c\\
# y(x) = \exp \left[-\int p(x)dx \right]\left(\int \exp \left[\int p(x) dx \right] r(x)dx + c \right)\\
# \end{align}
# \begin{align}
# y(x) = e^{-\int p dx} \left[\int e^{\int p dx} r dx + c \right] && \rightarrow \text{general solution for linear $1^\circ$ ODE}
# \end{align}
# 

# * If homogeneous, then $r=0$ and $y(x) = ce^{-\int p dx}$ as we found before.
# 

# ### Example
# $\frac{dy}{dx} = (x+1)^2 - y$

# 
# \begin{align}
# y' + y = (x+1)^2 && \rightarrow \text{non-homogeneous linear $1^\circ$ ODE}
# \end{align}
# $\therefore$ we know we can solve by turning it into an exact equation using an I.F. of the form 
# \begin{align}
# F(x) = \exp\left(\int p(x)dx \right) && \text{where $p(x) = 1$ (the coeff. of y)}\\
# F(x) = \exp\left(\int 1 dx \right) = e^x
# \end{align}
# * Multiply our ODE by our I.F.:\
# \begin{align}
# e^xy' + e^x y &= e^x(x+1)^2\\
# \frac{d}{dx}(e^xy) &= e^x(x+1)^2\\
# e^xy &= \int e^x (x+1)^2 dx\\
# e^xy &= e^x (x^2+1) + c && \text{integration by parts, often worst step. See this calc in next section}\\
# y(x) &= (x^2+1) + ce^{-x} && \text{general solution}
# \end{align}
# * Or we could have found this from our fully generalized solution
# \begin{align}
# y(x) &= e^{-\int p dx} \left[\int e^{\int p dx} r dx + c\right], && \text{where $p(x) = 1$ and $r(x) = (x+1)^2$}\\
# y(x) &= e^{-x}\left[\int e^x (x+1)^2 dx + c \right]\\
# &= (x^2 + 1) + ce^{-x}
# \end{align}

# # Integrating by Parts (Refresher)

# Method to integrate product of 2 terms
# \begin{align}
# \int u \cdot \frac{dv}{dx} dx = u \cdot v - \int v \cdot \frac{du}{dx} dx
# \end{align}
# 
# ### Example 
# Consider $\int e^x (x+1)^2 dx$
# 

# 
# $\rightarrow$ easiest to set 
# \begin{align}
# u = (x+1)^2 = x^2 + 2x + 1 && \text{and} && \frac{dv}{dx} = e^x
# \end{align}
# then,
# \begin{align}
# \frac{du}{dx} = 2x + 2 && \text{and} && v = e^x
# \end{align}
# 
# \begin{align}
# \int e^x (x+1)^2 dx = (x+1)^2 e^x - \int e^x (2x+2) dx && \text{need IBP again!}\\
# u = 2x + 2 \implies \frac{du}{dx} = 2\\
# \frac{dv}{dx} = e^x \implies v = e^x
# \end{align}
# 
# \begin{align}
# \text{Integral} &= (x+1)^2 e^x - \left[(2x+2)e^x - \int e^x \cdot 2 dx \right]\\
# &=  (x+1)^2 e^x - (2x+2)e^x + 2e^x\\
# &= e^x (x^2 + 2x + 1 - 2x -2 + 2)\\
# \int e^x (x+1)^2 dx &= e^x (x^2+1)
# \end{align}

# ## Recap of linear $1^\circ ODEs:$
# \begin{align}
# y' + p(x)y &= r(x)\\
# r(x) &= 0 \implies \text{homogeneous and separable}\\
# r(x) &\neq 0 \implies \text{non-homogeneous}\\
# \text{I.F.} = F(x) &= \exp\left[\int p(x) dx \right]
# \end{align}
# * Chain rule provides a simple solution:
# \begin{align}
# h = \int p(x) dx\\
# y(x) = e^{-h}\left[\int e^h r(x) dx + c \right] && \text{general solution}
# \end{align}
# Integral can be difficult to solve

# ## Chemical Engineering Example

# * Ethyl acetate decomposition
# \begin{equation}
# \ce{CH_3COOC_2H5 + H_2O -> CH_3COOH + C_2H_5OH}\\
# \text{(A)} \hspace{2cm} \text{$\rightarrow$first order kinetics}
# \end{equation}
# 
# 
# <img src="https://docs.google.com/drawings/d/e/2PACX-1vRmfzZPZ4UbOqbcO_9oZcw2odqzqUX9YR-v2VyA6_GkKJqU7qP9H1hfBwE9rmVPaNMtfo9jgyrEiub9/pub?w=316&h=314">
# CSTR : Continuous Stirred Tank Reactor

# * At steady state, each drop of feed is immediately converted to the composition of product stream.
# * What is concentration profile of A before steady state?\
# Mass balance on A:
# \begin{align}
# \text{Accum = In - Out + Gen - Consump}\\
# V\frac{dC}{dt} = FC_{in} - F C(t) - kV C(t)\\
# \frac{dC}{dt} = \frac{F}{V}C_{in} - \frac{1}{V}(F+kV)C
# \end{align}
# $\implies$ We can solve 2 ways:
# 1. Separate
# 2. Treat as linear equation

# ### 1. Separate
# \begin{align}
# \int \frac{dC}{\frac{FC_{in}}{V} - \frac{(F+kV)}{V}C} = \int dt
# \end{align}
# Useful integral:
# \begin{align}
# \int(ax + b)^{-1}dx = \frac{1}{a}\ln (ax + b) + c
# \end{align}
# Take $a=-\frac{F+kV}{V}$ and $b=\frac{FC_0}{V}$
# \begin{align}
# \frac{-V}{(F+kV)}\ln \left[ \frac{FC_{in}}{V} - \frac{(F+kV)}{V}C\right] &= t + \alpha\\
# \ln \left[ \frac{FC_{in}}{V} - \frac{(F+kV)}{V}C\right] &= -\frac{(F+kV)}{V}t + \frac{-(F+kV)}{V}\alpha
# \end{align}
# Take the last term as $\beta$:
# \begin{align}
# \frac{FC_{in}}{V} - \frac{(F+kV)}{V}C &= \exp \left[\frac{-(F+kV)t}{V} + \beta \right]\\
# &=\phi \exp \left[\frac{-(F+kV)}{V}t\right]\\
# \end{align}
# \begin{align}
# C(t) &= \frac{V}{(F+kV)}\left(\frac{FC_{in}}{V} - \phi \exp \left[\frac{-(F+kV)}{V} t \right] \right)\\
# C(t) &= \frac{FC_{in}}{F+kV} - \gamma \exp \left[\frac{-(F+kV)}{V} t \right] && \text{General Solution}
# \end{align}

# * Find $\gamma$ with I.C. $C(t=0) = C_0$ which is what the CSTR is charged with; doesn't have to be equal to $C_{in}$
# \begin{align}
# C(0) = C_0 &= \frac{FC_{in}}{F+kV} - \gamma \exp(0)\\
# \gamma &= \frac{FC_{in}}{F+kV} - C_0
# \end{align}
# \begin{align}
# C(t) &= \frac{FC_{in}}{F+kV} - \left(\frac{FC_{in}}{F+kV}-C_0 \right) \exp \left[-\frac{(F+kV)}{V}t \right]\\
# C(t) &= \frac{FC_{in}}{F+kV}\left(1 - \exp \left[-\frac{(F+kV)}{V}t \right] \right) + C_0 \exp \left[-\frac{(F+kV)}{V}t \right] && \text{specific  solution}
# \end{align}
# 
# Check : Does $C(0) = C_o$? (yes)\
# Does $C(t)$ as $T\rightarrow \infty$ = $\frac{F}{(F+kV)}C_{in}$? (yes!)

# In[ ]:




