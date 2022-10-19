#!/usr/bin/env python
# coding: utf-8

# {{ badge }}
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$

# # Recap for coupled differential equations:
# Coupled homogeneous 1st order ODEs:
# * $\vec{y'} = \arr{A} \vec{y}$ where $\vec{y}$ is a vector of unknown functions. 
# * $\vec{y} = [y_1, y_2, y_3, ... y_n]^T$ represents an $n^{th}$ order system with $n\times n \ \arr{A}$ 
# * We assumed that $\vec{x}e^{\lambda t}$ was a solution, then showed that it was. 
# * We proved that the eigenvalues and eigenvectors of $\arr{A}$ yield the solution
# \begin{align}
# \vec{y} = c_1\vec{x}^{(1)}e^{\lambda_1t} + c_2\vec{x}^{(2)}e^{\lambda_2t} + ... + c_n\vec{x}^{(n)}e^{\lambda_n t}
# \end{align}
# * In the last class, we looked at examples of coupled $1^\circ$ ODEs with real, distict eigenvalues.
# 

# ## Now let's look at complex eigenvalues
# 
# $\underline{\text{Ex}}$:
# \begin{align}
# y_1' &= 2y_2\\
# y_2' &= -2y_1\\
# \begin{bmatrix}y_1' \\ y_2' \end{bmatrix} &= \begin{bmatrix} 0 & 2 \\ -2 & 0 \end{bmatrix} \begin{bmatrix} y_1 \\ y_2 \end{bmatrix}\\
# \implies \vec{y}' &= \arr{A}\vec{y} \text{  where } \arr{A} = \begin{bmatrix} 0 & 2 \\ -2 & 0 \end{bmatrix}
# \end{align}
# 
# Can you guess what the solutions are going to look like? 

# 1. Find eigenvalues of $\arr{A}$:
# \begin{align}
# \det(\arr{A} - \lambda\arr{I}) = 0\\
# \left| \begin{array}{}0 - \lambda & 2 \\ -2 & 0-\lambda \end{array}\right| = 0 = \lambda^2 + 4\\
# \lambda = \pm 2i \rightarrow \lambda_1 = 2i, \  \lambda_2 = -2i
# \end{align}
# 2. Find eigenvectors\
# Set $(\arr{A} - \lambda\arr{I})\vec{x}^{(1)} = \vec{0}$
#   * For $\lambda_1 = 2i$
#   \begin{align}
#   \begin{bmatrix} -2i & 2 \\ -2 & -2i\end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}\\
#   \implies -2ix_1 + 2x_2 = 0\\
#   [-2x_1 - 2ix_2 = 0]i && \text{(redundant)}
#   \rightarrow -2ix_1 + 2x_2 = 0\\
#   ix_1 = x_2 \implies \vec{x}^{(1)} = \begin{bmatrix} 1 \\ i \end{bmatrix}
#   \end{align}
#   * For $\lambda_2 = -2i$ 
#   \begin{align}
# \begin{bmatrix} 2i & 2 \\ -2 & 2i \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}\\
# \implies 2ix_1 + 2x_2 = 0\\
# -ix_1 = x_2 \implies \vec{x}^{(2)} = \begin{bmatrix} 1 \\ i \end{bmatrix} = \begin{bmatrix} i \\ 1 \end{bmatrix}
#   \end{align}
# $\therefore$ General solution is:
# \begin{align}
# \vec{y} = c_1 \begin{bmatrix} 1 \\ i \end{bmatrix}e^{2it} + c_2 \begin{bmatrix} 1 \\ -i \end{bmatrix} e^{-2it}
# \end{align}
# Not very useful since it's imaginary
# * We will again use Euler's Formula:
# \begin{align}
# e^{2it} &= \cos(2t) + i\sin(2t)\\
# e^{-2it} &= \cos(-2t) + i\sin(-2t)\\
# &= \cos(2t) - i\sin(2t)
# \end{align}
# * Looking at part of the solution:
# \begin{align}
# \begin{bmatrix} 1 \\ i \end{bmatrix} e^{2it} &= \begin{bmatrix} 1 \\ i \end{bmatrix} (\cos 2t + i\sin 2t)\\
# &= \begin{bmatrix} \cos 2t + i\sin 2t \\ i\cos 2t - \sin 2t \end{bmatrix}\\
# &= \begin{bmatrix} \cos 2t \\ -\sin 2t \end{bmatrix} + i \begin{bmatrix} \sin 2t \\ \cos 2t\end{bmatrix}\\
# &= \vec{a} + i\vec{b}\\
# \begin{bmatrix} 1 \\ -i \end{bmatrix} e^{-2it} &= \begin{bmatrix} 1 \\ -i \end{bmatrix} (\cos 2t - i\sin 2t)\\
# &= \begin{bmatrix} \cos 2t - i\sin 2t \\ -i\cos 2t - \sin 2t \end{bmatrix}\\
# &= \begin{bmatrix} \cos 2t \\ -\sin 2t \end{bmatrix} - i \begin{bmatrix} \sin 2t \\ \cos 2t\end{bmatrix}\\
# &= \vec{a} - i\vec{b}
# \end{align}

# * Can we form real solutions from the imaginary solutions?\
# $\implies$ Let's try multiplying sum of imaginary solutions by $\frac{1}{2}$
# \begin{align}
# \frac{1}{2} \left\{\begin{bmatrix} 1 \\ i \end{bmatrix} e^{2it} + \begin{bmatrix} 1 \\ -i\end{bmatrix}e^{-2it} \right\} &= \frac{1}{2}(\vec{a} + i\vec{b} + \vec{a} - i\vec{b})\\
# &= \frac{1}{2}(2\vec{a})\\
# &= \vec{a} && \rightarrow\text{the sum of two linearly independent solutions is also a solution}\\
# -\frac{i}{2} \left\{\begin{bmatrix} 1 \\ i \end{bmatrix} e^{2it} - \begin{bmatrix} 1 \\ -i\end{bmatrix}e^{-2it} \right\} &= -\frac{i}{2}(\vec{a} + i\vec{b} - \vec{a} + i\vec{b})\\
# &= -\frac{i}{2}\cdot 2i\vec{b}\\
# &= \vec{b} && \rightarrow\text{the difference of two linearly independent solutions is also a solution}
# \end{align}
# $\therefore \vec{y} = c_1 \vec{a} + c_2 \vec{b}$ is a real representation of the general solution
# \begin{align}
# \vec{y} = c_1 \begin{bmatrix} \cos 2t \\ -\sin 2t \end{bmatrix} + c_2 \begin{bmatrix} \sin 2t \\ \cos 2t \end{bmatrix}
# \end{align}
# OR 
# \begin{align}
# y_1(t) &= c_1 \cos 2t + c_2 \sin 2t\\
# y_2(t) &= -c_1\sin 2t + c_2\cos 2t
# \end{align}
# 
# ## General solution for complex eigenvalues
# 
# We don't have to go through all of this work every time! 
# 
# 
# * For pairs of complex conjugate eigenvalues $\lambda\pm  \omega i$ with eigenvectors $\vec{x}=\vec{a}\pm \vec{b}i$ (the eigenvectors are also complex conjugates), the form will be:
# \begin{align}
# \vec{y}=e^\lambda t\left[c_1(\vec{a}\cos\omega t-\vec{b}\sin \omega t)+c_2(\vec{a}\sin\omega t+\vec{b}\cos \omega t)\right]
# \end{align}
# * To be clear, here $\vec{a}$ is the real part of the eigenvector, and $\vec{b}$ is the imaginary part of the eigenvector. 
# 
# * For the above example ($\lambda=\pm 2i$, eigenvector $\vec{x}=[1,\pm i]$, this yields
# \begin{align}
# \vec{y}&=e^0 t\left[c_1\left(\begin{bmatrix}1\\0\end{bmatrix}\cos2 t-\begin{bmatrix}0\\1\end{bmatrix}\sin 2 t\right)+c_2\left(\begin{bmatrix}1\\0\end{bmatrix}\sin2 t+\begin{bmatrix}0\\1\end{bmatrix}\cos 2 t\right)\right]\\
# y_1&=c_1\cos2t +c_2 \sin2t\\
# y_2&=-c_1\sin 2t+c_2\cos2t
# \end{align}

# ### Numerical solution to this example
# 
# Remember, when we solve this numerically, we also need initial conditions. So let's say that $y_1(t=0)=2$ and $y_2(t=0)=0.5$

# In[1]:


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

A = np.array([[0, 2],
              [-2,0]])

def linear_homogeneous_diffeq(t, y):
  return A@y

def linear_diffeq(t, y):
  return [2*y[1],
          -2*y[0]]


tspan = [0, 10]
y0 = [2, 0.5]
t_eval = np.linspace(0,10,100)

sol = solve_ivp(linear_homogeneous_diffeq, tspan, y0, t_eval=t_eval)

plt.plot(sol.t, sol.y.T,'o-')
plt.xlabel('Time t')
plt.ylabel('y')


# ## More than two coupled linear equations:
# 
# Very easy to add more to our systems
# 
# 
# $\underline{\text{Ex}}$: 
# \begin{align}
# y_1' &= 3y_1 + 5y_2 + 3y_3\\
# y_2' &= 4y_2 + 6y_3\\
# y_3' &= y_3
# \end{align}
# * Note that equation 3 isn't coupled, so in this case we can either solve that separately and plug into first two, or treat as 3D system:
# \begin{align}
# \vec{y}' = \arr{A}\vec{y} \ \text{ where } \vec{y} = \begin{bmatrix} y_1 \\ y_2 \\ y_3 \end{bmatrix} \text{ and } \arr{A} = \begin{bmatrix} 3 & 5 & 3 \\ 0 & 4 & 6 \\ 0 & 0 & 1 \end{bmatrix}
# \end{align}
# * The solution will be of the form:
# \begin{align}
# \vec{y} = c_1 \vec{x}^{(1)}e^{\lambda_1 t} + c_2 \vec{x}^{(2)}e^{\lambda_2 t} + c_3 \vec{x}^{(3)}e^{\lambda_3 t} 
# \end{align}
# where $\lambda$'s are the eigenvalues of $\arr{A}$ and $\vec{x}^{(i)}$'s are the eigenvectors of $\arr{A}$
# \begin{align}
# \arr{A} &= \begin{bmatrix} 3 & 5 & 3 \\ 0 & 4 & 6 \\ 0 & 0 & 1 \end{bmatrix} \\
# \implies \lambda_1 &= 3,\ \vec{x}^{(1)} = \begin{bmatrix} 1 & 0 & 0 \end{bmatrix}^T\\
# \lambda_2 &= 4,\ \vec{x}^{(2)} = \begin{bmatrix} 5 & 1 & 0 \end{bmatrix}^T\\
# \lambda_3 &= 1,\ \vec{x}^{(3)} = \begin{bmatrix} 7 & -4 & 2 \end{bmatrix}^T
# \end{align}
# \begin{align}
# \therefore \vec{y} = c_1 \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix} e^{3t} + c_2 \begin{bmatrix} 5 \\ 1 \\ 0 \end{bmatrix} e^{4t} + c_3 \begin{bmatrix} 7 \\ -4 \\ 2 \end{bmatrix} e^{t} 
# \end{align}
# OR:
# \begin{align}
# y_1(t) &= c_1 e^{3t} + 5c_2e^{4t} + 7c_3e^t\\
# y_2(t) &= c_2e^{4t} - 4c_3e^t\\
# y_3(t) &= 2c_3e^t 
# \end{align}
# * It's easy to check the solution to $y_3$ as it can be decoupled and is a solution to $y_3' = y_3$. 
# * Will need three initial conditions to find particular solution.
# 

# In[2]:


import numpy as np

# Quick numerical check to convince us that the eigenvalues/eigenvectors 
# are correct above
A = np.array([[3,5,3],
              [0,4,6],
              [0,0,1]])

eigval,eigvec = np.linalg.eig(A)

print(eigval)
print(eigvec[:,2])
print(eigvec[:,2]*7/eigvec[0,2])


# ### Numerical solution
# We need three initial conditions now. Let's say $y_1(0)=4, y_2(0)=2, y_3(0)=6$. Almost no changes from the above code

# In[3]:


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

A = np.array([[3, 5, 3],
              [0, 4, 6],
              [0, 0, 1]])

def linear_homogeneous_diffeq(t, y):
  return A@y

tspan = [0, 1]
y0 = [4,2,6]
t_eval = np.linspace(0,1,100)

sol = solve_ivp(linear_homogeneous_diffeq, tspan, y0, t_eval=t_eval)

plt.plot(sol.t, sol.y.T,'o-')
plt.xlabel('Time t')
plt.ylabel('y')


# 
# ## How do we handle non-homogeneous systems 
# 
# $\rightarrow$ Same idea as single equations
# \begin{align}
# y_1' &= -3y_1 + y_2 + 3\cos t\\
# y_2' &= y_1 -3y_2 -2\cos t -3\sin t\\
# \begin{bmatrix} y_1' \\ y_2' \end{bmatrix} &= \begin{bmatrix} -3 & 1 \\ 1 & -3 \end{bmatrix} \begin{bmatrix} y_1 \\ y_2 \end{bmatrix} + \begin{bmatrix} 3\cos t \\ -2\cos t - 3\sin t \end{bmatrix}\\
# \vec{y}' &= \arr{A}\vec{y} + \vec{g}(t) \text{ where $\vec{g}$ is a vector with no dependence on $\vec{y}$}
# \end{align}
# 
# 1. Solve homogeneous system. From $\arr{A}$:
# \begin{align}
# \lambda_1 = -4 && \vec{x}^{(1)} = \begin{bmatrix} -1 & 1 \end{bmatrix}^T\\
# \lambda_2 = -2 && \vec{x}^{(2)} = \begin{bmatrix} 1 & 1 \end{bmatrix}^T
# \end{align}
# \begin{align}
# \therefore \vec{y}_H (t) = c_1\begin{bmatrix} -1 \\ 1 \end{bmatrix} e^{-4t} + c_2 \begin{bmatrix} 1 \\ 1 \end{bmatrix} e^{-2t}
# \end{align}
# 2. Solve non-homogeneous part
#   * Now, the particular solution $\vec{y}_P$ is a vector
#   * For this example use MoUC. Choose $\vec{y}_P$ based on $\vec{g} = \begin{bmatrix} 3\cos t \\ -2\cos t -3\sin t \end{bmatrix}$
#   \begin{align}
#   \vec{y}_P = \vec{u}\sin t + \vec{v}\cos t
#   \end{align}
#   where $\vec{u} = \begin{bmatrix} u_1 \\ u_2 \end{bmatrix}$ and $\vec{v} = \begin{bmatrix} v_1 \\ v_2 \end{bmatrix}$ are the vectors of undetermined coefficients.\
#   Then, $\vec{y}_P' = \vec{u} \cos t - \vec{v}\sin t$
#   * Plug in expressions for $\vec{y}_P$ and $\vec{y}_P'$ into $\vec{y}_P' = \arr{A}\vec{y}_P + \vec{g}$
#   \begin{align}
#   \vec{u}\cos t - \vec{v}\sin t = \arr{A}(\vec{u}\sin t + \vec{v}\cos t) + \vec{g}\\
#   (\vec{u} - \arr{A}\vec{v})\cos t - (\arr{A}\vec{u} + \vec{v})\sin t = \vec{g}
#   \end{align}
#   * Match coefficients:
#     * cos(t) terms (two equations):
#   \begin{align}
#     \vec{u} - \arr{A}\vec{v} = \begin{bmatrix} 3 \\ -2 \end{bmatrix} 
#     \end{align}
#     * sin(t) terms (two equations):
#     \begin{align}
#     -\arr{A} \vec{u} -\vec{v} = \begin{bmatrix} 0 \\ -3 \end{bmatrix} 
#     \end{align}
#   * Solve system:\
#   From first equation: $\vec{u} = \begin{bmatrix} 3 \\ -2 \end{bmatrix} + \arr{A}\vec{v}$\
#   Plug into second equation:
#   \begin{align}
#   -\arr{A}\left(\begin{bmatrix} 3 \\ -2 \end{bmatrix} + \arr{A}\vec{v} \right) - \vec{v} = \begin{bmatrix} 0 \\ -3 \end{bmatrix}\\
#   \text{simplifies to }\ \begin{bmatrix} -11 & 6 \\ 6 & -11 \end{bmatrix} \begin{bmatrix}v_1 \\ v_2\end{bmatrix} = \begin{bmatrix} -11 \\ 6 \end{bmatrix}\\
#   \left[\begin{array}{ll|l}-11 & 6 & -11 \\ 6 & -11 & 6\end{array}\right] \text{ G.E} \rightarrow \left[\begin{array}{ll|l} 1&0&1 \\ 0&1&0\end{array} \right]\\
#   \therefore \vec{v} = \begin{bmatrix} 1 \\ 0 \end{bmatrix} \text{ and } \vec{u} = \begin{bmatrix} 3 \\ -2 \end{bmatrix} + \arr{A}\begin{bmatrix} 1 \\ 0 \end{bmatrix} = \begin{bmatrix} 0 \\ -1 \end{bmatrix}\\
#   \therefore \vec{y}_P = \begin{bmatrix} 0 \\ -1 \end{bmatrix} \sin t + \begin{bmatrix} 1 \\ 0 \end{bmatrix} \cos t
#   \end{align}
#   and
# \begin{align}
# \vec{y}(t) = c_1 \begin{bmatrix} -1 \\ 1 \end{bmatrix} e^{-4t} + c_2 \begin{bmatrix} 1 \\ 1 \end{bmatrix} e^{-2t} + \begin{bmatrix} 0 \\ -1 \end{bmatrix} \sin t + \begin{bmatrix} 1 \\ 0 \end{bmatrix} \cos t
# \end{align}
# 
# Notice that this says the initial conditions don't matter after a while!! 

# ### Numerical solutions

# In[4]:


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

A = np.array([[-3, 1],
              [1,-3]])

def linear_nonhomogeneous_diffeq(t, y):
  return A@y + [3*np.cos(t),
                -2*np.cos(t)-3*np.sin(t)]

tspan = [0, 20]
y0 = [-100, 400000]
t_eval = np.linspace(0,20,100)

sol = solve_ivp(linear_nonhomogeneous_diffeq, tspan, y0, t_eval=t_eval)

plt.plot(sol.t, sol.y.T,'o-')
plt.xlabel('Time t')
plt.ylabel('y')
plt.ylim([-2,2])


# ## Non-homogeneous system - VOP Method
# * Find the general solution to $\vec{y}' = \begin{bmatrix}-2 & 1 \\ 1 & -2  \end{bmatrix}\vec{y} + \begin{bmatrix} 2e^{-t} \\ 3t \end{bmatrix}$\
# same as:
# \begin{align}
# y_1' &= -2y_1 + y_2 + 2e^{-t}\\
# y_2' &= y_1 -2y_2 + 3t\\
# \vec{y}' &= \arr{A}\vec{y} + \vec{g}(t)
# \end{align}
# 1. Solve homogeneous system\
# from $\arr{A} \rightarrow |\arr{A}-\lambda\arr{I}|=0$
# \begin{align}
# \left|\begin{array}{} -2-\lambda & 1 \\ 1 & -2-\lambda \end{array}\right| = 0 & = (-2-\lambda)(-2-\lambda) -1\\
# &= +4 + 2\lambda + 2 \lambda + \lambda^2 - 1\\
# &= \lambda^2 + 4\lambda + 3 = (\lambda + 3)(\lambda + 1)\\
# \lambda_1 &= -3 ; \ \lambda_2 = -1
# \end{align}
# then $(\arr{A} - \lambda\arr{I})\vec{x} = \vec{0}$\
# For $\lambda_1=-3$ 
# \begin{align}
# \left[\begin{array}{ll|l} 1 & 1 & 0 \\ 1 & 1 & 0 \end{array} \right] \rightarrow x_1 + x_2 = 0 \implies x_1 = -x_2 \implies \vec{x}^{(1)} = \begin{bmatrix} 1 \\ -1 \end{bmatrix}
# \end{align}
# For $\lambda_2=-1$ 
# \begin{align}
# \left[\begin{array}{ll|l} -1 & 1 & 0 \\ 1 & -1 & 0 \end{array} \right] \rightarrow -x_1 + x_2 = 0 \implies x_1 = x_2 \implies \vec{x}^{(2)} = \begin{bmatrix} 1 \\ 1 \end{bmatrix}
# \end{align}
# Then,
# \begin{align}
# \vec{y}_H(t) = c_1 \begin{bmatrix} 1 \\ -1 \end{bmatrix} e^{-3t} + c_2\begin{bmatrix} 1 \\ 1 \end{bmatrix}e^{-t}
# \end{align}

# 2. Solve non-homogeneous part using Variation of Parameters
# * Write the homogeneous solution as:
# \begin{align}
# \vec{y}_H(t) &= \begin{bmatrix} y_1^{(1)} & y_1^{(2)} \\ y_2^{(1)} & y_2^{(2)} \end{bmatrix} \begin{bmatrix} c_1 \\ c_2 \end{bmatrix} = \arr{Y}(t)\vec{c}\  \text{ where $\arr{Y}(t)$ is the $``$fundamental matrix"}\\
# &= \begin{bmatrix} e^{-3t} & e^{-t} \\ -e^{-3t} & e^{-t}\end{bmatrix} \begin{bmatrix} c_1 \\ c_2 \end{bmatrix}
# \end{align}
# * Just like with one dimensional case, we will form a particular solution by "replacing" the arbitrary constants of the homogeneous solution with an unknown function of $t$, but this time with a vector. Assume $\vec{y}_P(t) = \arr{Y}(t)\vec{u}(t)$
# * We will determine $\vec{u}(t)$ by substituting $\vec{y}_P(t)$ into the original system \
# $\rightarrow$ let's derive in general terms

# \begin{align}
# \vec{y}' &= \arr{A} \vec{y} + \vec{g} \\
# \vec{y}_P' &= \arr{A}\vec{y}_P + \vec{g}\\
# (\arr{Y}\vec{u})' &= \arr{A} \arr{Y}\vec{u} + \vec{g}\\
# \arr{Y}'\vec{u} &+ \arr{Y}\vec{u}' = \arr{A}\arr{Y}\vec{u} + \vec{g}
# \end{align}
# $\implies$ because $\arr{Y}$ is a solution to the homogeneous system, we know it satifies $\arr{Y}' = \arr{A}\arr{Y}\rightarrow$ substitute this in 
# \begin{align}
# \arr{A}\arr{Y}\vec{u} +\arr{Y}\vec{u}' &= \arr{A}\arr{Y}\vec{u} + \vec{g}\\
# \arr{Y}^{-1}(\arr{Y}\vec{u}' &= \vec{g}) 
# \end{align}
# $\rightarrow$ can multiply with inverse. We know $\arr{Y}$ is non singular because $\det(\arr{Y})=$ Wronskian which is non-zero for a basis
# \begin{align}
# \vec{u}' &= \arr{Y}^{-1} \vec{g}\\
# \vec{u}(t) &= \int \arr{Y}^{-1}(t) \vec{g}(t)dt + \vec{c} && \rightarrow \text{each element is integrated separately}  
# \end{align}
# and 
# \begin{align}
# \vec{y}_P(t) &= \arr{Y}(t)\vec{u}(t)\\
# &= \arr{Y}(t) \int \arr{Y}^{-1}(t)\vec{g}(t) dt
# \end{align}
# and 
# \begin{align}
# \vec{y}(t) &= \vec{y}_H + \vec{y}_P\\
# \vec{y}(t) &= \arr{Y}(t)\vec{c} + \arr{Y}(t) \int \arr{Y}^{-1}(t)\vec{g}(t) dt
# \end{align}

# * Now for our example:\
# we need $\arr{Y}^{-1}(t)$. Can find using matrix property for a 2x2 matrix:
# \begin{align}
# \arr{A}^{-1} = \frac{1}{\det \arr{A}}\begin{bmatrix} a_{22} & -a_{12} \\ -a_{21} & a_{11} \end{bmatrix}
# \end{align}
# then 
# \begin{align}
# \arr{Y}^{-1} &= \frac{1}{e^{-4t} + e^{-4t}} \begin{bmatrix} e^{-t} & -e^{-t} \\ e^{-3t} & e^{-3t} \end{bmatrix}\\
# &=\frac{1}{2} \begin{bmatrix} e^{3t} & -e^{3t} \\ e^t & e^t \end{bmatrix}\\
# \arr{Y}^{-1}\vec{g} &= \frac{1}{2}\begin{bmatrix} e^{3t} & -e^{3t} \\ e^t & e^t \end{bmatrix} \begin{bmatrix} 2e^{-t} \\ 3t\end{bmatrix} \\
# &= \frac{1}{2} \begin{bmatrix} 2e^{2t} - 3te^{3t} \\ 2 + 3te^t \end{bmatrix}\\
# \vec{u}(t) = \int \arr{Y}^{-1}\vec{g} &= \frac{1}{2} \begin{bmatrix} \int 2e^{2t}dt - \int 3te^{3t}dt \\ \int 2 dt + \int 3t e^t dt \end{bmatrix}\\
# &=\frac{1}{2} \begin{bmatrix} e^{2t} - (te^{3t} - \frac{1}{3} e^{3t}) \\ 2t + 3(te^t - e^t) \end{bmatrix}
# \end{align}

# Then
# \begin{align}
# \vec{y}_P(t) &= \arr{Y}(t)\vec{u}(t)\\
# &= \begin{bmatrix} e^{-3t} & e^{-t} \\ -e^{-3t} & e^{-t} \end{bmatrix} \cdot \frac{1}{2} \cdot \begin{bmatrix} e^{2t}-te^{3t} + \frac{1}{3}e^{3t} \\ 2t + 3te^t - 3e^t \end{bmatrix} \\
# &= \frac{1}{2} \begin{bmatrix} e^{-t} -t + \frac{1}{3} + 2te^{-t}+ 3t -3 \\ -e^{-t} + t -\frac{1}{3} + 2te^{-t} + 3t - 3 \end{bmatrix}\\
# &= \begin{bmatrix} \frac{1}{2} e^{-t} + te^{-t} + t - \frac{4}{3} \\ -\frac{1}{2} e^{-t} + te^{-t} + 2t - \frac{5}{3}\end{bmatrix}\\
# \vec{y}_P(t) &= \frac{1}{2}\begin{bmatrix} 1 \\ -1 \end{bmatrix} e^{-t} + \begin{bmatrix} 1 \\ 1\end{bmatrix} te^{-t} + \begin{bmatrix} 1 \\ 2 \end{bmatrix} t -\frac{1}{3}\begin{bmatrix} 4 \\ 5 \end{bmatrix}
# \end{align}
# Then, 
# \begin{align}
# \vec{y}(t) &= \vec{y}_H(t) + \vec{y}_P(t)\\
# &= c_1 \begin{bmatrix} 1 \\ -1 \end{bmatrix} e^{-3t} + c_2 \begin{bmatrix} 1 \\ 1 \end{bmatrix} e^{-t} + \frac{1}{2}\begin{bmatrix} 1 \\ -1 \end{bmatrix} e^{-t} + \begin{bmatrix} 1 \\ 1 \end{bmatrix}te^{-t} + \begin{bmatrix} 1 \\ 2 \end{bmatrix} t - \frac{1}{3} \begin{bmatrix} 4 \\ 5 \end{bmatrix} && \text{general solution} 
# \end{align}

# ## Numerical solution
# 
# Integrate $\vec{y}' = \begin{bmatrix}-2 & 1 \\ 1 & -2  \end{bmatrix}\vec{y} + \begin{bmatrix} 2e^{-t} \\ 3t \end{bmatrix}$, with the initial condition $\vec{y}(t=1)=[4,3]^T$ from t=1 to t=10.

# In[5]:


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

A = np.array([[-2, 1],
              [1,-2]])

def linear_nonhomogeneous_diffeq(t, y):
  return A@y + [2*np.exp(-t),
                3*t]

tspan = [1, 10]
y0 = [4, 3]
t_eval = np.linspace(1,10,100)

sol = solve_ivp(linear_nonhomogeneous_diffeq, tspan, y0, t_eval=t_eval)

plt.plot(sol.t, sol.y.T,'o-')
plt.xlabel('Time t')
plt.ylabel('y')
# plt.ylim([-2,2])


# In[ ]:




