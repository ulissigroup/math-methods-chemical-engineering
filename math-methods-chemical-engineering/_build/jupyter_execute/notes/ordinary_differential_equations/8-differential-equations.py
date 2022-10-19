#!/usr/bin/env python
# coding: utf-8

# {{ badge }}
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$

# # Differential Equations
# 

# * See chapters 1-5 of text
# * A differential equation is an equation containing one or more detivatives of an unknown function.\
#   A) Ordinary Differential Equation (ODE)\
#     $\rightarrow$ Unknown depends on only one variable\
#     \begin{align}
#     \frac{dy}{dx} = x \rightarrow y=y(x)\\
#     \frac{dC}{dt} = -kC^2 \rightarrow C = C(t)
#     \end{align}
#   B) Partial Differential Equation (PDE)\
#   $\rightarrow$ Unknown function depends on two or more variables
#   \begin{align}
#   \frac{\partial y}{\partial t} = k \frac{\partial y}{\partial x} \rightarrow y=y(t,x)
#   \end{align}

# ## Notation
# * For an unknown function y(x)
#   1. First derivative of $y$: $\frac{dy}{dx} \equiv y'$
#   2. Second derivative: $\frac{d^2y}{dx^2}\equiv y''$
#   3. Third derivative: $\frac{d^3y}{dx^3}\equiv y'''$
# 
# ## Examples of ODEs:
# *   $y'=cos(x)$
# *   $y'' + 4y = 0$
# *   $x^2y'''y' + 2e^xy''=x^2+2$
# 
# "Order" of ODE: the order of the highest order derivative\
# $\rightarrow$ above equations are $1^{st}$, $2^{nd}$ and $3^{rd}$ order respectively.
# 

# 
# ## First Order Differential Equations
# \begin{align}
# F(x,y,y')=0 && \text{or} && y'=f(x, y)
# \end{align}
# will have a general solution
# \begin{align}
# y=h(x) && \rightarrow \text{explicit}\\
# H(x,y) = 0 && \rightarrow \text{implicit}
# \end{align}
# 
# ## Linear differential equations
# 
# A linear ODE is one where we can write it in the form $$a_n(t)y^{(n)}(t)+a_{n-1}(t)y^{(n-1)}(t)+\dots+a_1(t)y'(t)+a_0(t)y(t)=g(t)$$ where $y^n(t)$ is the nth derivative of y. If we can get it in this form, it is linear. If not, non-linear.
# 
# 

# ### Examples:
# * $y'=\cos(t)$
#   * Linear!
# * $y''+4y=0$ 
#   * Linear!
# * $x^2y'''y'+2e^xy''=x^2+2$
#   * Nonlinear because of the $y'''y'$ term
# * $\exp(y)+y''=4x$
#   * Nonlinear because of the $\exp(y)$ term
# * $\exp(x)+y''=3x$
#   * Linear
# * $y+y''=3x$
#   * Linear

# # Classification of Differential Equations
# 
# * We just discussed two types of differential equations:
#   * Partial Differential Equations (PDE's)
#   * Ordinary Differential Equations (ODE's)
# * We can also describe the order of the differential equation
#   * 1st order
#   * 2nd order 
#   * etc
# * Linearity of differential equation
#   * Linear
#   * Non-linear
# 
# As we progress in the class we will add more classifications.
# 

# ### Separable $1^\circ $ ODEs
# $\rightarrow$ should already know these\
# $\rightarrow$ Equation is separable if it can be written as:
# \begin{align}
# g(y)y'=f(x)
# \end{align}
# $\rightarrow$ since $y'=\frac{dy}{dx}$, then 
# \begin{align}
# g(y)dy = f(x)dx && \text{separate all $y$'s from all $x$'s}
# \end{align}
# Simple integration solution
# \begin{align}
# \int g(y)dy = \int f(x) dx + \underline{c} && \text{don't forget the constant!}
# \end{align}

# ### Example separable first order ODE 
# 
# $y'=2x$
# 

# \begin{align}
# \frac{dy}{dx} = 2x\\
# dy = 2x dx && \text{separate}\\
# \int dy =\int 2x dx + c && \text{integrate}\\
# y(x) = x^2 + c && \text{general solution}
# \end{align}
# We need more information to determine $c$. Once $c$ is known $\rightarrow$ particular solution.\
# e.g. if we know $y(x=0)=2$, then
# \begin{align}
# 2&=0^2+c\\
# c&=2
# \end{align}
# and the particular solution is $y(x) = x^2 + 2$
# * Separable $1^\circ$ ODE's arise frequantly in engineering problems

# ### Example separable first order ODE 
# 
# 
# $\underline{Ex}$: There is a pond just north of Pittsburgh where people like to swim in summer. Unfortunately on Memorial Day (May 25), the pond is polluted with a substance $T$ at a concentration of 100 ppm. There is a flow rate of fresh water into the pond that equals the rate at which water leaches into the soil (pond volume of $10^4 m^3$ doesn't change). Swimming will be safe only once $C_T < 5$ ppm. Will swimmers be back in business before Labor Day (Sept. 5)? 
# 
# <img src="https://docs.google.com/drawings/d/e/2PACX-1vS8mkzTjLckrqsayNi9kZhbf6jhU_D4Gw3g2Gz_x924UWb2MPNObOqqS7oqIl4j7Lpsc-yWrT6SbpYH/pub?w=960&h=720">
# 
# $\rightarrow$ Solve by performing a mass balance on what? (the pollutant)
# \begin{align}
# \text{Accumulation = In - Out + Gen - Consumption}
# \end{align}
# Inlet term is zero because only fresh water enters. Generation and Consumption are zero as well because there's no reaction.\
# $\rightarrow$ Unlike linear systems, this problem is interested in what happens before steady state.\
# $\rightarrow$ There is an element of time involved
# \begin{align}
# \frac{dM}{dt} = \frac{d(C \cdot V)}{dt} = -FC
# \end{align}
# Check units,
# \begin{align}
# \frac{kg}{s} [=] \frac{\frac{kg}{L} L}{s} [=] \frac{L}{s} \cdot \frac{kg}{L} 
# \end{align}

# $\rightarrow$ since $V$ is constant, we can pull it out
# \begin{align}
# V\frac{dC}{dt} = -FC && \text{$1^\circ$ ODE that models our system}\\
# \frac{dC}{C} = -\frac{F}{V} dt && \text{separate}\\
# \int \frac{dC}{C} = -\frac{F}{V} \int dt + \alpha && \text{integrate}\\
# \ln C = -\frac{F}{V} t + \alpha\\
# C(t) = \exp\left[-\frac{F}{V}t + \alpha\right] && \text{General solution}
# \end{align}
# $\rightarrow$ to find $\alpha$, we need to know $C(t)$ at some time $\rightarrow$ will yield our particular solution\
# $\rightarrow$ we know $C(t=0)=C_0=100ppm$
# \begin{align}
# \therefore C(t=0) = C_0 &= \exp\left[-\frac{F}{V} \cdot 0 + \alpha\right]\\
# C_0 &= \exp(\alpha)\\
# \alpha &= \ln C_0 && \rightarrow \text{can plug into general solution}\\
# \end{align}
# 
# \begin{align}
# C(t) &= \exp\left[-\frac{F}{V}t + lnC_0 \right]\\
# &=\exp\left[-\frac{F}{V}t\right]\cdot \exp[lnC_0]\\
# C(t) &= C_0 \exp \left[ -\frac{F}{V}t\right] && \text{particular solution}
# \end{align}

# For our problem: $C_0 =100ppm$, $V=10^4 m^3$, $F=100 m^3$/day
# \begin{align}
# C(t) &= 100 ppm \cdot \exp\left(\frac{100 m^3/day}{10^4m^3}t \right)\\
# C(t) &
# = 100 \exp(-0.01t) && \text{where $C[=]ppm, t[=]days$} 
# \end{align}
# $\rightarrow$ Can swimmers swim on Labor Day? Let's determine how long it will take for $C(t)$ to drop to 5 ppm.
# \begin{align}
# 5=100\exp(-0.01t)\\
# ln\left(\frac{5}{100} \right) = -0.01t\\
# t = 300 days
# && \rightarrow \text{Can't swim until next summer :(}
# \end{align}

# # Transient Batch Reactor 
# <img src="https://docs.google.com/drawings/d/e/2PACX-1vSf2jVdqbl_bfEDKJRTToTinA-WikP4MJbZv8E-t8GKjxe9CNSQMRdO9Vinb2YBDz8qLcpk1Fi5I2iq/pub?w=960&h=720">
# \begin{equation}
# \ce{CH_3COOC_2H5 + H_2O -> CH_3COOH + C_2H_5OH}\\
# \text{ethyl acetate}\hspace{3cm} \text{acetic acid} \hspace{0.8cm} \text{ethanol}
# \end{equation}
# 
# 
# 

# $\rightarrow$ Develop a model for the concentration of ethyl acetate in the reactor as a function of time. The reaction consumes E.A. at a rate proportional to the instantaneous concentration of E.A.\
# $\rightarrow$ Assume the reactor is charged with $C_0$.\
# Step 1: What is changing with time?\
# Concentration of E.A. $\rightarrow$ mass is conserved.\
# MB on E.A: Accum = In - Out + Gen - Consump\
# The reactor has no inlet and outlet. E.A. isn't being produced by the reaction. Hence, Generation term is also zero. 
# \begin{align}
# \frac{dM}{dt} = -kC
# \end{align}
# We were told Consump $\propto C_{EA}$. Given $C[=]\frac{mass}{vol} \implies k[=]\frac{vol}{time}$ to match LHS units of $\frac{mass}{time}$
# 
# 

# \begin{align}
# V\frac{dC}{dt} = -kC\\
# \frac{dC}{C} = -\frac{k}{V} dt && \text{separate}\\
# \ln C = -\frac{k}{V}t + \beta\\
# C(t) = \exp\left(-\frac{k}{V}t \right) \cdot \exp(\beta)
# \end{align}
# Let $\exp(\beta) = \alpha$,
# \begin{align}
# C(t) = \alpha \cdot \exp\left(-\frac{k}{V}t\right) && \text{general solution}
# \end{align}
# $\rightarrow$ We know that $C(t=0) = C_0$
# \begin{align}
# \therefore C(t=0) = C_0 = \alpha \cdot \exp(0) \implies C_0 = \alpha
# \end{align}
# then, 
# \begin{align}
# C(t) = C_0 \cdot \exp\left(-\frac{k}{V} t\right) && \text{particular solution}
# \end{align}
# 
# * What if, instead of knowing the initial concentration, we were told that ater 10 minutes, conc. = $C_1$?\
# $\rightarrow$ only $\alpha$ will change
# \begin{align}
# C(t=10) = C_1 = \alpha \cdot \exp\left(\frac{-10k}{V} \right) \implies \alpha=\frac{C_1}{\exp\left(\frac{-10k}{V} \right)}
# \end{align}
# and
# \begin{align}
# C(t) = \frac{ C_1 \exp\left(-\frac{k}{V} t \right)}{\exp\left(-\frac{10k}{V} t \right)} =C_1 \exp\left[\frac{k}{V} (10 - t)\right] 
# \end{align}

# The basic batch reactor pops up often.\
# $\rightarrow$ Redo the problem with different kinetics:\
# What if E.A. is consumed at a rate proportional to the square of the instantaneous concentration of E.A.?\
# Now: Accumulation = - Consumption
# \begin{align}
# V\frac{dC}{dt} = -kC^2 && \text{where now $k[=]\frac{vol^2}{mass \cdot time}$}\\
# \frac{dC}{C^2} = -\frac{k}{V}dt\\
# -\frac{1}{C} = -\frac{k}{V}t + \alpha\\
# \frac{1}{C} = \frac{k}{V}t + \alpha && \text{OK to change sign of $\alpha$ since it is arbitrary}\\
# C(t) = \frac{1}{\frac{k}{V}t+ \alpha} && \text{general solution}
# \end{align}
# 
# Given $C(0) = C_0$, then $C(t=0) = 1/\alpha \implies \alpha = \frac{1}{C_0}$
# \begin{align}
# C(t) = \frac{1}{\frac{k}{V}t+ \frac{1}{C_0}}\\
# C(t) = \frac{C_0}{\frac{k}{V}C_0t+ 1} && \text{particular solution}
# \end{align}

# # Reality checks
# 
# * Consider the possible steady states of the system. 
# * What do you know? What will be a reasonable answer?
# * If you're worried about your solution, plug it in and check if it's a solution

# In[ ]:




