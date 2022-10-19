#!/usr/bin/env python
# coding: utf-8

# {{ badge }}

# # The Shooting Method for Solving BVPs
# 
# Let's go back to that pirate defense problem from last class. We figured out that it wasn't possible to hit the pirate ship with a cannonball speed of 40 m/s, regardless of how we aimed the cannon. Let's try and solve this again using the same method we used before for finding steady states to a problem, `root`.
# 
# \begin{align}
# \frac{d}{dt}\begin{bmatrix}x\\z\\v_x\\v_z \end{bmatrix} = 
# \begin{bmatrix} v_x\\v_z\\-rv_x\\-rv_z-g\end{bmatrix} \rightarrow \text{physics with air friction}
# \end{align}
# \begin{align}
# x(t=0) = 0\\
# z(t=0) = H\\
# z(x=L) = 0 
# \end{align}
# 
# 

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

g = 9.8 #m/s^2
r = 0.1 #1/s

def diff_eq(t, y):
  x, z, vx, vz = y

  return [vx,
          vz,
          -r*vx,
          -r*vz-g]

tspan = [0, 10] # t=0s to 10s

v = 40 
theta = 5/360*2*np.pi
y0 = [0,100,v*np.cos(theta), v*np.sin(theta)]
t_eval = np.linspace(0,10)

sol = solve_ivp(diff_eq, tspan, y0, t_eval=t_eval)

# plot z vs x
plt.plot(sol.y[0,:], sol.y[1,:])
plt.plot(200,0,'ok')
plt.xlabel('Distance [m]')
plt.ylabel('Height [m]')




# ## Reformulate the ODE
# 
# As posed, this is a little difficult to solve, since we don't know exactly how long to integrate for to hit the pirate ship. Let's solve for all of the variables in terms of derivatives with respect to the horizontal position. 
# 
# \begin{align}
# \frac{d}{dx}\begin{bmatrix}x\\z\\v_x\\v_z \end{bmatrix} =\frac{dt}{dx}\frac{d}{dt}\begin{bmatrix}x\\z\\v_x\\v_z \end{bmatrix}= \frac{1}{v_x}\frac{d}{dt}\begin{bmatrix}x\\z\\v_x\\v_z \end{bmatrix}=\frac{1}{v_x}\begin{bmatrix} v_x\\v_z\\-rv_x\\-rv_z-g\end{bmatrix} =
# \begin{bmatrix} 1\\v_z/v_x\\-r\\(-rv_z-g)/v_x\end{bmatrix} 
# \end{align}
# 
# We don't even need the first equation anymore (of course $d/dx(x)=1$), so we're left with 
# \begin{align}
# \frac{d}{dx}\begin{bmatrix}z\\v_x\\v_z \end{bmatrix} =
# \begin{bmatrix} v_z/v_x\\-r\\(-rv_z-g)/v_x\end{bmatrix} 
# \end{align}
# 

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

g = 9.8 #m/s^2
r = 0.1 #1/s

def diffeq(x, y):
  z, vx, vz = y
  return [vz/vx,
          -r,
          (-r*vz-g)/vx]


v = 47.5 
theta = 45/360*2*np.pi
y0 = [100,v*np.cos(theta), v*np.sin(theta)]

sol = solve_ivp(diffeq,
                [0,200],
                y0,
                t_eval=np.linspace(0,200))

plt.plot(sol.t, sol.y[0,:])
plt.xlabel('Distance from shore [m]')
plt.ylabel('Height [m]')
plt.xlim([0,200])
plt.ylim([0,200])




# Great! Now all we have to do is figure out the final height at x=200, and make that equal to zero. First, let's just see what the final height is at x=200.

# In[ ]:


print('The final height is %1.2f m'%sol.y[0,-1]) 


# Let's see what happens when we try to use 40 m/s now. Try going even lower. What do you see? *What* happens?

# In[ ]:





# Great! Now all we have to do is figure out the final height at x=200, and make that equal to zero. First, let's just see what the final height is at x=200.

# In[ ]:





# ## Root finding to solve the BVP
# 
# Ok, now we're in a great position to find the exact angle or initial conditions that will solve this BVP. Let's remind ourselves how `fsolve` works. 
# 
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fsolve.html
# 
# fsolve will take a function $\vec{y}=f(\vec{x})$, and try to find the $\vec{x}$ such that $\vec{y}=0$. If we can phrase our problem in terms of this one, `fsolve` will give us the answer!
# 
# 

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

g = 9.8 #m/s^2
r = 0.1 #1/s

def diffeq(x, y):
  z, vx, vz = y
  return [vz/vx,
          -r,
          (-r*vz-g)/vx]


def BVP_satisfy(x):
  v, theta_deg = x
  theta_rad = theta_deg/360*2*np.pi
  y0 = [100,v*np.cos(theta_rad), v*np.sin(theta_rad)]
  
  sol = solve_ivp(diffeq,
                [0,200],
                y0,
                t_eval=np.linspace(0,200))

  print(f'{v:1.17f}, {theta_deg:1.17f}')
  return [sol.y[0,-1],
          0]

fsolve_sol = fsolve(BVP_satisfy,[30,50])

print(fsolve_sol)


# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

g = 9.8 #m/s^2
r = 0.1 #1/s

def diffeq(x, y):
  z, vx, vz = y
  return [vz/vx,
          -r,
          (-r*vz-g)/vx]


v = 30
theta = 50/360*2*np.pi
y0 = [100,v*np.cos(theta), v*np.sin(theta)]

sol = solve_ivp(diffeq,
                [0,200],
                y0,
                t_eval=np.linspace(0,200))

plt.plot(sol.t, sol.y[0,:])
plt.xlabel('Distance from shore [m]')
plt.ylabel('Height [m]')
plt.xlim([0,200])
plt.ylim([0,200])

sol.y[0,-1]


# We guessed 49 m/s, and it found that actually 48.6 m/s worked best for hitting the pirate ship. Let's plot the guess solution and the final solution. 

# In[ ]:





# There are many initial conditions that will satisfy the final boundary conditions in this problem. 

# In[ ]:





# ## Summary of the shooting method to solve BVPs
# 
# This method of solving BVPs is called the shooting method, because you guess initial conditions and shoot over to other values to check whether they work or now. The key steps are:
# 1. You have to be able to solve the IVP first, assuming you know the initial conditions. If you have mistakes in the IVP, you won't get any farther than that.
# 2. Once you have the IVP integration working, you can try adjusting the initial conditions a little to see what happens and make sure everything is working correctly.
# 3. When you are confident about your solution, write a function for `root` that takes in initial conditions you have to guess, and returns all zeros when the boundary conditions are satisfied. This could be integrating to a fixed point, and calculating how far the boundary conditions are from the desired values.
# 4. Use `root` to vary the initial conditions and find a solution. 
# 5. Try a few different guesses for the initial conditions, to see if there are multiple solutions available.
# 
# This is a more complex problem then we've done in the past, because it requires three different levels of numerics:
# * numpy/arrays/etc for coding the system of differential equations
# * `solve_ivp` to integrate those differential equations
# * `root` to adjust initial conditions and solve the problem.
# 
# There are other methods of solving BVPs numerically which use a grid of points from the initial to final independent variable. We won't cover those here, but you can look into them with `solve_bvp`. It will handle the boundary condition solution for you, but the downside is that you have to worry about a mesh of points. 
# 
# 

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

g = 9.8 #m/s^2
r = 0.1 #1/s

def diffeq(x, y):
  z, vx, vz = y
  return [vz/vx,
          -r,
          (-r*vz-g)/vx]


def BVP_satisfy(x):
  v, theta_deg = x
  theta_rad = theta_deg/360*2*np.pi
  y0 = [100,v*np.cos(theta_rad), v*np.sin(theta_rad)]
  
  sol = solve_ivp(diffeq,
                [0,200],
                y0,
                t_eval=np.linspace(0,200),
                atol=1e-8,
                rtol=1e-8)

  # print(f'{v:1.17f}, {theta_deg:1.17f}')
  return [sol.y[0,-1],
          0]

final_params = [fsolve(BVP_satisfy,[30,theta_guess]) for theta_guess in np.linspace(-20,50)]

stacked_parameters = np.stack(final_params)

plt.plot(stacked_parameters[:,1],stacked_parameters[:,0],'o')


# # Chemical Engineering Example
# 

# <img src="https://docs.google.com/drawings/d/e/2PACX-1vSf2jVdqbl_bfEDKJRTToTinA-WikP4MJbZv8E-t8GKjxe9CNSQMRdO9Vinb2YBDz8qLcpk1Fi5I2iq/pub?w=106&h=182">
# \begin{align}
# 2A \leftrightarrow B\\
# B \rightarrow 2C
# \end{align}
# At $t=0$, operator add 5 mol of $A$ into a 5L vessel that contained an unknown amount of $B$ and $C$ beforehand. After 5 minutes, we stopped the reactor and measured
# \begin{align}
# C_B = 0.02 M\\
# C_C = 2 M
# \end{align}
# 
# How much $B$/$C$ was initially in tank?\
# How much A left at the end?
# 
# First, let's get this into a set of differential equations
# \begin{align}
# \frac{dC_A}{dt} &= -2k_1C_A^2+2k_{-1}C_B\\
# \frac{dC_B}{dt} &= k_1C_A^2-k_{-1}C_B - k_2C_B\\
# \frac{dC_C}{dt} &= 2k_2C_B 
# \end{align}
# 
# We also need the boundary conditions.
# \begin{align}
# C_A(t=0) &= 1M\\
# C_B(t=t_f) &= 0.02M\\
# C_C(t=t_f) &= 2 M\\
# k_1 =1 [1/M/min], k_2 &= 1[1/min], k_{-1}=0.2[1/min]
# \end{align}
# 
# This is a set of nonlinear coupled ODE's, so it would be hard to solve analytically. Our best best is to use numerical methods, but the boundary conditions are specified at both t=0 and $t=t_f$. We're going to have to use `solve_ivp` and `root`.
# 

# In[ ]:





# Now, specify a function that will take in one vector $x$ containing the initial conditions for $C_B$ and $C_C$ and return how far the final conditions for those concentrations are from the known final values. Use this with `root` to find the initial conditions and solve for the final concentration profile. 

# In[ ]:





# In[ ]:





# In[ ]:




