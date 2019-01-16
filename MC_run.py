#!/usr/bin/env python
# coding: utf-8

# ## Test of accelerated Monte Carlo ##

# In[1]:


from particle import particle
from physics import em_physics
from cylinder import cylinder

import matplotlib.pyplot as plt
import numpy as np
from numba import jit
import pandas as pd

import json

#get_ipython().run_line_magic('matplotlib', 'inline')


# ### Define the Geometry ###

# In[2]:


# define the geometry
radius = 50 # cm
height = 100 #cm
cryostat = cylinder(R=radius,h=height)

radius = 40 # cm
height = 80 #cm
fiducial = cylinder(R=radius,h=height)


# ### Setup the photon physics ###

# In[3]:


# define the physics
em = em_physics()


# ### Event generation ###

# In[ ]:


#
# number of events
#
nevent = 1000001
#
# gamma energy
#
energy = 1000 # keV
#
# output filename
#
mcout = 'mcdata.json'

#
# open output data file
#
f = open(mcout,'w')
f.write('[\n')

#
# event generation loop
#
for ieve in range(nevent):
    if ieve%25000 == 0:
        print("generated ",ieve," events")
    #
    # make a particle
    #
    p = particle(type='gamma',
                 energy=energy, 
                 geometry=cryostat, 
                 fiducial=fiducial, 
                 vrt='fiducial_scatter',
                 edep_max=250,
                 nscatter_max=1,
                 physics=em,
                 debug=False)
    
    #
    # propagate the particle and retrieve the intersection points
    #
    p.propagate()

    #
    # store data
    #
    for i in range(len(p.xint)):
        #
        # data structure for each energy deposit
        #
        dd = {'x':p.xint[i][0][0],'y':p.xint[i][0][1],'z':p.xint[i][0][2],'de':p.xint[i][1],
              'w':p.weight,'n':p.nscatter}
        json.dump(dd,f, sort_keys=True, indent=4)
        f.write(',\n')
         
#
# close output data file
#
f.write('{}]')
f.close()


# In[ ]:





# In[ ]:




