#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/ulissigroup/math-methods-chemical-engineering/blob/master/lecture_notes/07b-numerical-eigenvalues.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
# $$\newcommand{\arr}[1]{\underline{\underline{#1}}}$$    
# $$\newcommand{\vec}[1]{\underline{#1}}$$   
# $$\require{mhchem}$$

# # Eigenvectors in python
# 
# Unsurprisingly, there is a function to calculate eigenvalues and eigenvectors in python! 
# * For most cases, we can use the ```np.linalg.eig``` function
# * If we only wanted the eigenvalues, ```np.linalg.eigvals``` will just calculate those. 
# 

# ## Example 1: real roots
# 
# Let's start with a simple example where we already know the roots: a lower triangular matrix
# \begin{align}
# \arr{A}=\begin{bmatrix}6&10&6 \\ 0&8&12 \\ 0&0&2 \end{bmatrix} 
# \end{align}

# In[1]:


import numpy as np

A = np.array([[6,10,6],
              [0,8,2],
              [0,0,2]])

eigenvalues, eigenvectors = np.linalg.eig(A)

print(eigenvalues)
print(eigenvectors)


# In[ ]:
































# In[2]:


import numpy as np

A = np.array([[6,10,6], 
              [0,8,12],
              [0,0,2]])

# calculate the eigenvalues of A:
eigenvalues = np.linalg.eigvals(A)
print('The eigenvalues of A are: %s'%eigenvalues)

#Calculate both the eigenvalues and eigenvectors of A
eigenvalues, eigenvectors = np.linalg.eig(A)
# calculate the eigenvalues of A:
print('The eigenvalues of A are: %s'%eigenvalues)

#Print the first eigenvalue and eigenvector
print('The eigenvalue corresponding to lambda=%d is: %s'%(eigenvalues[0],eigenvectors[:,0]))

# We can also iterate over each eigenvalue/eigenvector pair
for i in range(len(eigenvalues)):
  print('The eigenvalue corresponding to lambda=%d is: %s'%(eigenvalues[i],eigenvectors[:,i]))


# Let's compare these to the known values from lecture 6:
# 
# Eigenspace corresponding to $\arr{A}$ :
# \begin{array}{}
# 6, \begin{bmatrix} 1\\0\\0 \end{bmatrix} ; & 8,\begin{bmatrix} 5\\1\\0 \end{bmatrix}; & 2, \begin{bmatrix} 7\\-4\\2\end{bmatrix} ; & \vec{x}=\vec{0}
# \end{array}
# 
# Notice that the eigenvectors returned by numpy are the same ratios, but different absolute numbers. We can rescale these if we need. 
# 

# In[3]:


eigenvectors[:,1]/eigenvectors[1,1]


# In[ ]:





































# In[4]:


print('The eigenvalue corresponding to lambda=%d is: %s'%(eigenvalues[2],eigenvectors[:,2]/eigenvectors[2,2]*2))


# ## Complex roots
# 
# $ \arr{A} = \begin{bmatrix} 1&2 \\ -2&1 \end{bmatrix} $
# 
# This is easy now that we know how to use ```np.linalg.eig```!

# In[5]:


A = np.array([[1,2],
              [-2,1]])

eigenvalues, eigenvectors = np.linalg.eig(A)

print(eigenvalues)
print(eigenvectors)


# In[6]:


(1+2.j)**2


# In[ ]:
































# In[7]:


A = np.array([[1,2],[-2,1]])
eigenvalues, eigenvectors = np.linalg.eig(A)

print(eigenvalues)
print(eigenvectors)


# Notice that numpy is using "j" to indicate a complex number. Otherwise things look pretty simple!
# 
# ### Here are some examples of complex numbers:

# In[ ]:































# In[8]:


a = 1+2.j
print('a=%s'%str(a))
print('a+1=%s'%str(a+1))
print('2a=%s'%str(2*a))


# # In-class problem 
# 
# Calculate the eigenvalues and eigenvectors for the matrix 
# 
# $\arr{A} = \begin{bmatrix} 0&9&-12 \\ -9&0&20 \\ 12&-20&0 \end{bmatrix}$
# 

# # Population eigenvalue problem
# 
# Let's start with calculating the eigenvalues and eigenvectors of the population transition matrix $\arr{A}$

# In[ ]:





# In[ ]:






















# In[9]:


A = np.array([[0.3,0.2,0],
              [0.6,0.2,0],
              [0.1,0.6,1]])

eigenvalues, eigenvectors = np.linalg.eig(A)

print(eigenvalues)
print(eigenvectors)


# Matches up with what we already know! 
# 
# Let's do a little example for a random starting population of 1000 people. 

# In[ ]:





# In[ ]:



























# In[10]:


# Generate three random numbers
p = np.random.rand(3)

#Normalize the random numbers and multiple by 1000 so that they all add correctly.
p = np.round(p/p.sum()*1000)

print('Starting population: %d alive, %d sick, %d dead '%(p[0],p[1],p[2]))


# Now, let's increment the population:

# In[ ]:





# In[11]:


p2 = A@p
print('Starting population: %d alive, %d sick, %d dead '%(p2[0],p2[1],p2[2]))


# In[12]:


p2 = np.linalg.matrix_power(A,20)@p
print('Starting population: %d alive, %d sick, %d dead '%(p2[0],p2[1],p2[2]))


# Let's do this after n iterations:
# 

# In[ ]:





# In[ ]:

































# In[13]:


p2 = np.linalg.matrix_power(A,20)@p
print('Starting population: %d alive, %d sick, %d dead '%(p2[0],p2[1],p2[2]))


# Notice that everyone dies. This happens for any starting configuration. 

# ## In-class problem: zombies
# 
# We discover that there's a 20% chance of a dead person coming back from the dead (a zombie). Calculate the steady state distribution of live, sick, and dead people. 

# In[14]:


A = np.array([[0.3,0.2,0.0],
              [0.6,0.2,0.2],
              [0.1,0.6,0.8]])

eigenvalues, eigenvectors = np.linalg.eig(A)

print(eigenvalues)
print(eigenvectors)


# In[15]:


p2 = np.linalg.matrix_power(A,20)@p
print('Starting population: %d alive, %d sick, %d dead '%(p2[0],p2[1],p2[2]))


# In[ ]:





# In[ ]:




