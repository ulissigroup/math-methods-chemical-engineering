#!/usr/bin/env python
# coding: utf-8

# {{ badge }}
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$
# 

# # Recap on PDEs from last lecture
# 
# * PDE's are more complicated than ODE's
# * You need to analyze each independent variable (x, y, t, etc), determine the order, and then determine the number of boundary / initial conditions you need
# * Many PDE's that are important in chemical engineering are linear and second order
#   * $u_t=c^2u_{xx}$
#   * $u_{tt}=c^2u_{xx}$
#   * $\nabla^2u=0$
# 
# * Separation of variables for linear PDEs:
#   1. Assume $u(x,t)=X(x)T(t)$
#   2. Plug in to PDE, separate, set equal to a constant to yield multiple ODE's
#   3. Solve each ODE
#   4. Get general solution
#   5. Apply boundary conditions to constrain the problem
#   6. Initial conditions require a linear combination of the possible solutions with coefficients $c_n$ that we talked about last time. 
#   

# ## Reminder: constant temperature initial condition for heated bar
# 
# 
# Solving a problem with a bar that has a sin temperature profile is not very interesting. How would you make a bar like that? I have no idea.
# 
# Let's do a more reasonable set of initial conditions $f(x)=T_0=5$, that is, we immerse the bar in a water bath 5 K warmer than our hands, then take it out and hold it out on both sides. We don't have time to show how to solve this equation right now (we'll do that next class). The magic $c_n$ that specify this are 
# \begin{align*}
# c_n&=
# \begin{cases} 
#       \frac{4T_0}{n\pi}  & n\text{ odd} \\
#       0 & n\text{ even}
#          \end{cases}
# \end{align*}
# That makes our final solution
# \begin{align*}
# T(x,t)&=\sum_{n=1}^{\infty} c_ne^{-\lambda_n t} \sin\left(\frac{n\pi x}{L}\right)\\
# T(x,t)&=\sum_{n=1}^{\infty} \frac{20L}{n\pi}e^{-\lambda_n t} \sin\left(\frac{n\pi x}{L}\right) \text{ if n odd}
# \end{align*}
# 
# with 
# \begin{align*}
# \lambda_n=\alpha\left(\frac{n\pi}{L}\right)^2
# \end{align*}

# In[1]:


import numpy as np
import matplotlib.pyplot as plt

from matplotlib import animation, rc
from IPython.display import HTML

L=1 #m
alpha=0.1 #m^2/s
lambdan=alpha*(np.pi/L)**2
T0 = 5

x = np.linspace(0, L, 200)

# First set up the figure, the axis, and the plot element we want to animate
fig, ax = plt.subplots()
plt.close()

T0 = 5*np.ones(x.shape)
ax.plot(x,T0,'--',label='Initial Temp Profile')

ax.set_xlim(( -0.5, 1.5))
ax.set_ylim((0, 6))
ax.set_xlabel('x')
ax.set_ylabel('Temperature [K]')
line, = ax.plot([], [], lw=2, label='Temperature Profile')
ax.legend()

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return (line,)

# animation function. This is called sequentially  
def animate(i):
  T = np.zeros(len(x))
  for n in range(1,100):
    if n%2==1:
      cn=4*T0/np.pi/n
    else:
      cn=0
    lambdan=alpha*(n*np.pi/L)**2
    T=T+cn*np.exp(-lambdan*t[i])*np.sin(np.pi*n*x/L)
    
  line.set_data(x, T)
  return (line,)
  
t = np.linspace(0,3,100)

anim = animation.FuncAnimation(fig, animate, init_func=init,
                             frames=100, interval=100, blit=True)

# Note: below is the part which makes it work on Colab
rc('animation', html='jshtml')
anim
  


# In summary, the initial value problem we're faced with is:
# \begin{align*}
# T(x,t) = f(x) = T_0 = \sum_{n=1}^\infty c_n e^{-\lambda_n 0} \sin\left(\frac{n\pi x}{L}\right)= \sum_{n=1}^\infty c_n  \sin\left(\frac{n\pi x}{L}\right)
# \end{align*}

# ## Fourier Series
# 
# The Fourier series is a special series approximation to a function in terms of sin and cos terms
# \begin{align*}
# f(x) = a_0 + \sum_{n=1}^\infty (a_n \cos nx + b_n \sin nx)
# \end{align*}
# where $f(x)$ is an arbitrary function and $c_n$ are the **Fourier coefficients**. 
# 
# * Because the Fourier series is in terms of sin/cos, it is periodic
# * This process is guaranteed to converge. 
# 
# 
# How do we find the Fourier coefficients $c_n$?
# * Special formula for sin/cos for each term:
# \begin{align*}
# a_0 &= \frac{1}{2\pi}\int_{-\pi}^{\pi}f(x)dx\\
# a_n &= \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\cos nx dx\\
# b_n &= \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\sin nx dx
# \end{align*}
# Notice that the range is defined from $-\pi$ to $\pi$ - we will deal with that later. 
# 
# ### Example: flat line $f(x)=c$
# 
# We need to evaluate each of those terms above. The constant term is easy
# \begin{align*}
# a_0 &= \frac{1}{2\pi}\int_{-\pi}^{\pi}f(x)dx\\
# &=\frac{1}{2\pi}\int_{-\pi}^{\pi}cdx=c
# \end{align*}
# Now let's do the sin/cos terms
# \begin{align*}
# a_n &= \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\cos nx dx\\
# &=\frac{1}{\pi}\int_{-\pi}^{\pi}c\cos nx dx=0\\
# b_n &= \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\sin nx dx\\
# &=\frac{1}{\pi}\int_{-\pi}^{\pi}c\sin nx dx=0
# \end{align*}
# So the Fourier series for $f(x)=c$ is just $f(x)=c$. Easy! It was harder for the bar example because we were contrained to $T(0,t)=T(L,t)=0$. 

# 
# ### Example: square  wave
# 
# Let's consider 
# \begin{align*}
# f(x)=\begin{cases}-k& -\pi < x<0\\
# k& 0<x<\pi
# \end{cases}
# \end{align*}
# 
# We still need to evaluate each of those terms above.
# \begin{align*}
# a_0 &= \frac{1}{2\pi}\int_{-\pi}^{\pi}f(x)dx\\
# a_0 &= \frac{1}{2\pi}\int_{-\pi}^{0}f(x)dx+\frac{1}{2\pi}\int_{0}^{\pi}f(x)dx\\
# a_0 &= \frac{1}{2\pi}\int_{-\pi}^{0}-kdx+\frac{1}{2\pi}\int_{0}^{\pi}kdx=0\\
# \end{align*}
# Now for the $a_n$ terms
# \begin{align*}
# a_n &= \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\cos nx dx\\
# &= \frac{1}{\pi}\int_{-\pi}^{0}f(x)\cos nx dx+\frac{1}{\pi}\int_{0}^{\pi}f(x)\cos nx dx\\
# &= \frac{1}{\pi}\int_{-\pi}^{0}-k\cos nx dx+\frac{1}{\pi}\int_{0}^{\pi}k\cos nx dx\\
# &= \frac{1}{\pi}\left[ \left. \frac{-k\sin nx}{n}\right|_{-\pi}^{0}+\left. \frac{k\sin nx}{n}\right|_{0}^{\pi}\right]\\
# &= \frac{1}{\pi}\left[ 0-0+0-0\right]=0\\
# \end{align*}
# So all the $a_n$ terms are zero. One last one.
# \begin{align*}
# b_n &= \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\sin nx dx\\
# &= \frac{1}{\pi}\int_{-\pi}^{0}f(x)\sin nx dx+\frac{1}{\pi}\int_{0}^{\pi}f(x)\sin nx dx\\
# &= \frac{1}{\pi}\int_{-\pi}^{0}-k\sin nx dx+\frac{1}{\pi}\int_{0}^{\pi}k\sin nx dx\\
# &= \frac{1}{\pi}\left[ \left. \frac{k\cos nx}{n}\right|_{-\pi}^{0}+\left. \frac{-k\cos nx}{n}\right|_{0}^{\pi}\right]
# \end{align*}
# This is a little tricky. There are two possibilities for n odd and n even. If n is odd
# \begin{align*}
# b_n= \frac{1}{\pi}\left[\frac{k}{n}-\frac{-k}{n}+\frac{k}{n}-\frac{-k}{n}\right]=\frac{4k}{\pi n}\\
# \end{align*}
# If n is even, then $\cos 0=\cos n\pi$ and $b_n=0$.
# So the Fourier series for this function is 
# \begin{align*}
# f(x) &= a_0 + \sum_{n=1}^\infty (a_n \cos nx + b_n \sin nx)\\
# &=\frac{4k}{\pi}\sin x + \frac{4k}{3\pi}\sin 3x + \dots
# \end{align*}
# Look familiar? Let's plot this!

# In[2]:


import numpy as np
import matplotlib.pyplot as plt

from matplotlib import animation, rc
from IPython.display import HTML

k=5

x = np.linspace(-2*np.pi , 2*np.pi, 200)

# First set up the figure, the axis, and the plot element we want to animate
fig, ax = plt.subplots()
plt.close()

f = -5*(x<0)+5*(x>0)
ax.plot(x,f,'--',label='f(x)')


ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.set_ylim((-k-2,k+2))

line, = ax.plot([], [], lw=2, label='Approximation to f(x)')
last_term, = ax.plot([], [], lw=2, label='Latest term')

ax.legend()

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return (line,)

# animation function. This is called sequentially  
def animate(i):
  f_approx = np.zeros(len(x))
  term = np.zeros(len(x))
  for n in range(1,n_max[i]+1):
    if n%2==1:
      cn=4*k/np.pi/n
    else:
      cn=0

    term = cn*np.sin(n*x)
    f_approx+=cn*np.sin(n*x)
    
  line.set_data(x, f_approx)
  last_term.set_data(x,term)
  ax.set_title('n=%d'%n_max[i])
  return (line,)
  

n_max = list(range(1,60))

anim = animation.FuncAnimation(fig, animate, init_func=init,
                             frames=len(n_max), interval=100, blit=True)

# Note: below is the part which makes it work on Colab
rc('animation', html='jshtml')
anim
  


# Notice the little wiggles around the corners of the square wave? This effect is called Gibbs phenomenon, which occur when a Fourier series is trying to accomodate a discontuinity. 
# 

# ## Special properties of the Fourier series
# 
# This idea doesn't work for any arbitrary functions. The sin/cos series has some special properties. The most important property is that the series is orthogonal. The integral of any two of the terms multiplied together is zero. 
# \begin{align*}
# \int_{-\pi}^\pi \cos nx \cos mx dx=0, \text{ if }n\neq m\\
# \int_{-\pi}^\pi \sin nx \sin mx dx=0, \text{ if }n\neq m\\
# \int_{-\pi}^\pi \sin nx \cos mx dx=0
# \end{align*}
# These are relatively simple to prove with a little bit of trigonometry. 
# 
# ## Other orthogonal series 
# 
# * Most sequences wouldn't have this property. However, there are other orthogonal series covered in section 11.6 of your textbook. These orthogonal series come up in engineering when solving various PDE's in different geometries (rectangular, cylindrical, spherical etc). 
# 
# * You might see these ideas in heat/mass transfer (if solving for heat transfer in a spherical particle), reaction engineering (reactions in a spherical particle), or quantum mechanics in physical chemistry. 
# 
# * Note that the form of the solution is coming from the ODE we solved in separation of variables, so this should match the basis functions you're using.
# 
# ## Tricks to make solutions easier 
# 
# Notice how we just did everything in terms of $-\pi$ to $\pi$ but the original PDE was from 0 to L? There are a few tricks to switch from period $2\pi$ to period 2L, and to simplify the fourier series. These are covered in section 11.2 in the textbook, and you'll have some simple examples on the next homework.
# 
# For example, if your function is odd ($f(x)=-f(-x)$), then you will get an expansion in terms of $\sin$. If your function is even $f(x)=f(-x)$, then you will get an expansion in terms of $\cos$. Notice that the bar example we have is actually an expansion for a square wave, that happens to work for $x=0$ to $x=\pi$ for the heated bar.

# ## Second example
# 
# Calculate $b_n$ in the Fourier series for $f(x)=|x|$ (the absolute value of x). You might find the following integrals helpful (from integration by parts):
# \begin{align*}
# \int x\sin axdx = \frac{\sin (ax)}{a^2}-\frac{x\cos ax}{a} && \int x\cos axdx = \frac{\cos (ax)}{a^2}+\frac{x\sin ax}{a}
# \end{align*}
# 
# * First, we solve for $a_0$
# 
# \begin{align*}
# a_0=\frac{1}{2\pi}\int_{-\pi}^\pi f(x)dx=\frac{1}{\pi}\int_0^\pi xdx=\frac{1}{\pi} \left[\frac{x^2}{2}\right]_0^\pi=\pi/2
# \end{align*}
# 
# * Next, we solve for $a_n$
# \begin{align*}
# a_n&=\int_{-\pi}^\pi f(x) \cos nx dx\\
# &=\int_{-\pi}^0 -x \cos nx dx+\int_{0}^\pi x \cos nx dx\\
# &=\frac{2}{\pi} \int_0^\pi x\cos nx dx\\
# &=\frac{2}{\pi}\left[\frac{\cos nx}{n^2}+\frac{x\sin nx}{n}\right]_0^\pi\\
# &=\frac{2}{\pi}\left[\frac{\cos n\pi}{n^2}-\frac{\cos 0}{n^2}\right]\\
# &=\frac{2(\cos n\pi-1)}{\pi n^2}
# \end{align*}
# 
# * Finally, $b_n$. This function is even, so right away we know that $b_n$ will be 0. However, you can still do the math and verify that. 
# 
# \begin{align*}
# b_n &= \frac{1}{\pi}\int_{-\pi}^{\pi}f(x)\sin nx dx\\
# &=\frac{1}{\pi}\int_{-\pi}^{0}-x\sin nx dx+\frac{1}{\pi}\int_{-\pi}^{\pi}x\sin nx dx\\
# &=\frac{1}{\pi}\left[-\left(\frac{\sin (nx)}{n^2}-\frac{x\cos nx}{n}\right)_{-\pi}^0 +\left(\frac{\sin (nx)}{n^2}-\frac{x\cos nx}{n}\right)_{0}^\pi \right]\\
# &=\frac{1}{\pi}\left[-\left(\frac{\sin (0x)}{n^2}-\frac{0\cos 0x}{n}\right)+\left(\frac{\sin (-n\pi)}{n^2}-\frac{-\pi\cos -n\pi}{n}\right)+\left(\frac{\sin (n\pi)}{n^2}-\frac{\pi\cos n\pi}{n}\right)-\left(\frac{\sin (n0)}{n^2}-\frac{0\cos n0}{n}\right)\right]\\
# &=\frac{1}{\pi}\left[-\left(\frac{\sin (0x)}{n^2}-\frac{0\cos 0x}{n}\right)+\left(\frac{\sin (-n\pi)}{n^2}-\frac{-\pi\cos -n\pi}{n}\right)+\left(\frac{\sin (n\pi)}{n^2}-\frac{\pi\cos n\pi}{n}\right)-\left(\frac{\sin (n0)}{n^2}-\frac{0\cos n0}{n}\right)\right]\\
# &=\frac{1}{\pi}\left[\frac{\pi\cos n\pi}{n}-\frac{\pi\cos n\pi}{n}\right]=0
# \end{align*}
# 
# We can also plot this series solution to verify it.
# 

# In[3]:


import numpy as np
import matplotlib.pyplot as plt

from matplotlib import animation, rc
from IPython.display import HTML


x = np.linspace(-2*np.pi , 2*np.pi, 200)

# First set up the figure, the axis, and the plot element we want to animate
fig, ax = plt.subplots()
plt.close()

f = np.abs(x)
ax.plot(x,f,'--',label='f(x)')

ax.set_xlabel('x')
ax.set_ylabel('f(x)')
ax.set_ylim((-k-2,k+2))

line, = ax.plot([], [], lw=2, label='Approximation to f(x)')
last_term, = ax.plot([], [], lw=2, label='Latest term')

ax.legend()

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return (line,)

# animation function. This is called sequentially  
def animate(i):
  f_approx = np.zeros(len(x))

  #Add in the first term
  a0 = np.pi/2
  f_approx+=a0

  #Add in the current cos term
  term = np.zeros(len(x))
  for n in range(1,n_max[i]+1):
    an=2/np.pi*((np.cos(n*np.pi)-1)/n**2)
    term = an*np.cos(n*x)
    f_approx+=an*np.cos(n*x)
    
  #Plot the most recent term and the full approximation
  line.set_data(x, f_approx)
  last_term.set_data(x,term)
  ax.set_title('n=%d'%n_max[i])
  return (line,)
  

n_max = list(range(1,60))

anim = animation.FuncAnimation(fig, animate, init_func=init,
                             frames=len(n_max), interval=100, blit=True)

# Note: below is the part which makes it work on Colab
rc('animation', html='jshtml')
anim
  


# # Numerical solutions to PDE's
# 
# Numerical solutions to PDE's can be very complicated, especially in multiple dimensions or complicated geometries. For example, fluid flow in a reactor with multiple mixing regions:
# 
# ![alt text](https://umich.edu/~elements/5e/web_mod/radialeffects/femlab_doc/images/tank_mixer12.png)
# 
# The governing engineering/physics is usually straightforward. However, setting up the boundary conditions and solving is very complicated. The standard approach is to use a small mesh of points will little finite volumes/elements, discretize the PDE into a linear problem and solve that. This is sufficiently complicated that commercial packages like `openfoam`, `comsol`, or `ansys` (remember ansys hall?) are used. ansys generally has tools for each type of PDE you want to solve. `comsol` you can type in the PDE more generally and solve. These are taught in the graduate course 06-663 Analysis and Modeling of Transport Phenomena. 
# 

# # Note - repeated roots in coupled ODE's
# 
# One thing we brushed over when solving linear coupled ODE's is what happens when you have a repeated root. For example:
# \begin{align*}
# \vec{y}'=\arr{A}\vec{x}=\begin{bmatrix}
# 0&1\\-1&-2
# \end{bmatrix}\vec{x}
# \end{align*}
# The eigenvalues of this matrix are $\lambda=-1,-1$, which is a repeated root. There is a single eigenvalue $$\vec{x}^{(1)}=\begin{bmatrix}-1\\1\end{bmatrix}$$ We don't have enough information to construct the solution. The approach in 2nd order ODE's was to tack a $t$ onto the solution
# \begin{align*}
# \vec{y}=c_1\begin{bmatrix}-1\\1\end{bmatrix}e^{-t}+c_2\begin{bmatrix}-1\\1\end{bmatrix}te^{-t}
# \end{align*}
# **but that does not work here!** 
# 
# A second solution is 
# \begin{align*}
# y_2=\vec{x}^{(1)}te^{\lambda t}+\vec{\eta}e^{\lambda t}
# \end{align*}
# where $\vec{\eta}$ is a generalized eigenvector that solves the problem:
# \begin{align*}
# (\vec{A}-\lambda \arr{I})\vec{x}^{(2)}=\vec{x}^{(1)}
# \end{align*}
# 
# Plugging in $\lambda$ and $\vec{x}^{(1)}$, we get
# \begin{align*}
# \begin{bmatrix}1&1\\-1&-1\end{bmatrix}\vec{x}^{(2)}=\begin{bmatrix}-1\\1\end{bmatrix}
# \end{align*}
# This tells us that $x_1+x_2=-1$, so we can pick any eigenvector that satisfies this. Let's choose $\vec{x}^{(2)}=[0,-1]$. Our final solution will be
# \begin{align*}
# \vec{y}=c_1\begin{bmatrix}-1\\1\end{bmatrix}e^{-t}+c_2\begin{bmatrix}-1\\1\end{bmatrix}te^{-t}+c_2\begin{bmatrix}0\\-1\end{bmatrix}e^{-t}
# \end{align*}

# In[ ]:




