#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 08:10:19 2019

@author: alex
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 13:43:22 2019

@author: alex
"""
import numpy as np

# mesh and plotting imports
import pymesh
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.spatial import distance
#graph imports
import networkx as nx

    
graph_edge_example = np.array([[0,1],[1,2],[2,3],[3,0],[0,4],[1,4],[2,4],[3,4]])
graph_edge_example = ((0,1,{'weight':1} ),(1,2, {'weight' : 1}), 
                       (2,3, {'weight' : 1}),
                       (3,0, {'weight' : 1}),
                              (0,4, {'weight':np.sqrt(2)/2}),
                              (1,4, {'weight':np.sqrt(2)/2}),
                              (2,4, {'weight':np.sqrt(2)/2}),
                              (3,4, {'weight' : np.sqrt(2)/2}))
nodes = np.array([[0,0],[1,0],[1,1],[0,1],[1/2,1/2]])

# unweighted
   
graph_edge = np.array([[0,1],[1,2],[2,3],[3,0],[0,4],[1,4],[2,4],[3,4],[0,2]])
graph_edge = ((0,1,{'weight':1} ),(1,2, {'weight' : 1}), 
              (0,2,{'weight': 0.01}),
                       (2,3, {'weight' : 1}),
                       (3,0, {'weight' : 1}),
                              (0,4, {'weight':np.sqrt(2)/2}),
                              (1,4, {'weight':np.sqrt(2)/2}),
                              (2,4, {'weight':np.sqrt(2)/2}),
                              (3,4, {'weight' : np.sqrt(2)/2}))

# Dijkstra
G = nx.from_edgelist(graph_edge)
G_triangle = nx.from_edgelist(graph_edge_example)

length, path = nx.multi_source_dijkstra(G,{0})
length_t, path_t = nx.multi_source_dijkstra(G_triangle,{0})


#    convert lengths into array for triplot


dist_t = np.zeros(5)
for key, value in length.items():
    print(key,value)
    dist_t[key]= value
    
faces = [[0, 4, 1],[1,4,2], [2,4,3],[3,4,0]]
#plot
fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
tcf = ax1.plot(nodes[0,0], nodes[0,1] ,'r*', markersize=11)
    
tcf = ax1.tricontourf(nodes[:,0], nodes[:,1], faces, dist_t)
fig1.colorbar(tcf)
plt.show()

dist_triangle = np.zeros(5)
for key, value in length_t.items():
    print(key,value)
    dist_triangle[key]= value
    
#plot
fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
tcf = ax1.plot(nodes[0,0], nodes[0,1] ,'r*', markersize=11)
    
tcf = ax1.tricontourf(nodes[:,0], nodes[:,1], faces, dist_triangle)
fig1.colorbar(tcf)
plt.show()
    
    
    
    
