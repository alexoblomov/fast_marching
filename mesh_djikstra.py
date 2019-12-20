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

# edge weights: to represent distances in irregular mesh
def create_weighted_graph(graph):
    weighted_graph = []
    for edge in graph[1]:
       weight =({'weight' : 1/(distance.euclidean(graph[0][edge[0]], graph[0][edge[1]]) + 1)})
       weighted_graph.append((edge[0],edge[1], weight))
    
    return weighted_graph

def add_connectivity(graph):
    return 0
# load and visualize mesh
mesh = pymesh.load_mesh("poly_4.off");
plt.triplot(mesh.vertices[:,0], mesh.vertices[:,1], mesh.faces, 'ko-', lw = 0.5, alpha=0.5, ms = 0.7)

# load connectivity data
#connected_nodes = 

# convert mesh to graph 
graph = pymesh.mesh_to_graph(mesh)
# add weights to account for irregularity of mesh
weighted_graph = create_weighted_graph(graph)
# create networkx graph
G = nx.from_edgelist(weighted_graph)


# load sources
sources = np.loadtxt("poly_boundary_4.txt")
src = np.where(np.isin(graph[0][:,0:2], sources[:,0:2])[:,0])[0]

#
#graph_edge_example = np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
#graph_edge_example = ((0,1,{'weight':6} ),(0,2, {'weight' : np.sqrt(34)}), 
#                       (0,3, {'weight' : np.sqrt(13)}),
#                       (1,2, {'weight' : np.sqrt(34)}),
#                              (1,3, {'weight': np.sqrt(13)}),
#                              (2,3, {'weight' : 3}))
#nodes = np.array([[0,0],[6,0],[3,5],[3,2]])

## Dijkstra
#G_triangle = nx.from_edgelist(graph_edge_example)
#length_t, path_t = nx.multi_source_dijkstra(G_triangle,{3})

length, path = nx.multi_source_dijkstra(G,set(src))


#TEST

#convert lengths into array for triplot

#dist_triangle = np.zeros(4)
#for key, value in length_t.items():
#    dist_triangle[key]= value
#    
#faces = [[0, 1, 3],[0,2,3], [1,2,3]]
##plot
#fig1, ax1 = plt.subplots()
#ax1.set_aspect('equal')
#tcf = ax1.plot(nodes[3,0], nodes[3,1] ,'r*', markersize=11)
#    
#tcf = ax1.tricontourf(nodes[:,0], nodes[:,1], faces, dist_triangle)
#
#plt.show()




#####

#convert lengths into array for triplot
# =============================================================================
# dist_list=[]
# temp = []
# for key, value in length.items():
#     temp = [key,value]
#     dist_list.append(temp)
#     
# dist = np.array(dist_list)
# =============================================================================
dist= np.zeros(np.shape(graph[0])[0]) 
for key, value in length.items():
    dist[key]= value

#plot
fig1, ax1 = plt.subplots()
ax1.set_aspect('equal')
for i in np.arange(len(src)):
    tcf = ax1.plot(mesh.vertices[src[i]][0], mesh.vertices[src[i]][1] ,'r*', markersize=11)

tcf = ax1.tricontourf(mesh.vertices[:,0], mesh.vertices[:,1], mesh.faces, dist)

plt.show()
