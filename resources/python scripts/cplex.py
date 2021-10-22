#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
try:
    import cplex
except:
    if hasattr(sys, 'real_prefix'):
        #we are in a virtual env.
        get_ipython().system('pip install cplex')
    else:
        get_ipython().system('pip install --user cplex')
        
try:
    import docplex.mp
except:
    if hasattr(sys, 'real_prefix'):
        #we are in a virtual env.
        get_ipython().system('pip install docplex')
    else:
        get_ipython().system('pip install --user docplex')

# first import the Model class from docplex.mp
from docplex.mp.model import Model

# create one model instance, with a name
m = Model(name='check_popularity')

# Read input from file
import csv

edges = []
weights = []
variables = []
Dictionary = {}
dummy_Dictionary = {}
pointer = 0


with open('C:/Users/SAI/SMFQ_Graphmatching/resources/HR2LQ/input6_new_output_3.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for row in csv_reader:
        ver1 = row[0]
        ver2 = row[1]
        wt = row[2]
        is_matched = row[3]
        
        edges.append(ver1+'_'+ver2)
        weights.append(int(wt))
        
        # create variable for each edge
        variables.append(m.continuous_var(name="y"+str(pointer+1)))
        
        # constraint : variable for any edge should be greater than or equal to 0
        m.add_constraint(variables[pointer] >= 0)

        if "dummy" not in ver1:
            if ver1 not in Dictionary:
                Dictionary[ver1] = []
            Dictionary[ver1].append(pointer) 
        else:
            if ver1 not in dummy_Dictionary:
                dummy_Dictionary[ver1] = []
            dummy_Dictionary[ver1].append(pointer)
            
        if "dummy" not in ver2:
            if ver2 not in Dictionary:
                Dictionary[ver2] = []
            Dictionary[ver2].append(pointer) 
        else:
            if ver2 not in dummy_Dictionary:
                dummy_Dictionary[ver2] = []
            dummy_Dictionary[ver2].append(pointer)
            
        pointer = pointer + 1
# print(edges)
# print(weights)
# print(variables)
# print(Dictionary)

for key in Dictionary:
    m.add_constraint(m.sum(variables[i] for i in Dictionary[key]) == 1)
#     print(m.sum(variables[i] for i in Dictionary[key]) == 1)
#     print(key, '->', Dictionary[key])

for key in dummy_Dictionary:
    m.add_constraint(m.sum(variables[i] for i in dummy_Dictionary[key]) <= 1)


import numpy as np
m.maximize(m.sum(variables[i] * int(weights[i]) for i in range(len(variables))))
# m.maximize(np.dot(variables, weights))
# print(np.dot(variables, weights))
# m.get_objective_expr()

m.print_information()
s = m.solve()
m.print_solution()
# m.solve_details()

