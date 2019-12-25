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
       weight =({'weight' : distance.euclidean(graph[0][edge[0]], graph[0][edge[1]])})

       weighted_graph.append((edge[0],edge[1], weight))
    
    return weighted_graph

def add_connectivity(graph, weighted_graph, connect):
    nodes_to_connect = np.where(np.isin(graph[0][:,0:2], connect[:,0:2])[:,0])[0]
    weight = ({'weight' : 100000})
    row = 0
    num_rows = np.shape(connect)[0]
    while row < num_rows:
        weighted_graph.append((nodes_to_connect[row], nodes_to_connect[row+1], weight))
        row = row + 2
        
    return weighted_graph

for i in np.arange(4,16,2):
    mesh_name = "poly_" + str(i) + ".off"
    mesh_boundary = "poly_boundary_" + str(i) + ".txt"
    mesh_connectivity = "poly_connectivity_" + str(i) + ".txt"
    # load and visualize mesh
    mesh = pymesh.load_mesh(mesh_name);
#    plt.triplot(mesh.vertices[:,0], mesh.vertices[:,1], mesh.faces, 'ko-', lw = 0.5, alpha=0.5, ms = 0.7)
    
    # convert mesh to graph 
    graph = pymesh.mesh_to_graph(mesh)
    # add weights to account for irregularity of mesh
    weighted_graph = create_weighted_graph(graph)
    
    # add connectivity
    connect = np.loadtxt(mesh_connectivity)
    weighted_graph_with_connectivity = add_connectivity(graph, weighted_graph.copy(), connect)
    
    # create networkx graph
    G = nx.from_edgelist(weighted_graph)
    G_prime = nx.from_edgelist(weighted_graph_with_connectivity)
    
    # load sources
    sources = np.loadtxt(mesh_boundary)
    src = np.where(np.isin(graph[0][:,0:2], sources[:,0:2])[:,0])[0]    
    
    # run Dijkstra
    length, path = nx.multi_source_dijkstra(G,set(src))
    length_prime, path_prime = nx.multi_source_dijkstra(G_prime, set(src))
    
    #convert lengths into array for triplot
    dist= np.zeros(np.shape(graph[0])[0]) 
    for key, value in length.items():
        dist[key]= value
    
    #plot
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    for i in np.arange(len(src)):
        tcf = ax1.plot(mesh.vertices[src[i]][0], mesh.vertices[src[i]][1] ,'r*', markersize=11)
    tcf = ax1.tricontourf(mesh.vertices[:,0], mesh.vertices[:,1], mesh.faces, dist)
    plt.title('Dijkstra without connectivity')
    plt.show()
    
    dist_prime= np.zeros(np.shape(graph[0])[0]) 
    for key, value in length_prime.items():
        dist_prime[key]= value
        
    fig1, ax1 = plt.subplots()
    ax1.set_aspect('equal')
    for i in np.arange(len(src)):
        tcf = ax1.plot(mesh.vertices[src[i]][0], mesh.vertices[src[i]][1] ,'r*', markersize=11)
    tcf = ax1.tricontourf(mesh.vertices[:,0], mesh.vertices[:,1], mesh.faces, dist_prime)
    plt.title('Dijkstra with connectivity')
    plt.show()
