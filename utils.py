#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 07:48:03 2019

@author: alex

"""
from scipy.spatial import distance
import numpy as np
# edge weights: to represent distances in irregular mesh
def create_weighted_graph(graph):
    weighted_graph = []
    for edge in graph[1]:
       weight =({'weight' : 100* distance.euclidean(graph[0][edge[0]], graph[0][edge[1]])})
       weighted_graph.append((edge[0],edge[1], weight))
    
    return weighted_graph

def add_connectivity(nodes, weighted_graph, connect):
    nodes_to_connect = np.where(np.isin(nodes[:,0:2], connect[:,0:2])[:,0])[0]
    weight = ({'weight' : 0.001})
    row = 0
    num_rows = np.shape(connect)[0]
    while row < num_rows:
        weighted_graph.append((nodes_to_connect[row], nodes_to_connect[row+1], weight))
        row = row + 2
        
    return weighted_graph

def dijkstra(graph, initial, bound_nodes = [], voronoi = False):
    far = set(graph[0])
    accepted = {}
    path = {}
    considered = {}

    for j, init_node in enumerate(initial):
        accepted[init_node] = 0
        far.remove(init_node)
        for edge in graph[1][init_node]:
            if edge not in initial and edge not in considered.keys():
                considered[edge] = graph.distances[(init_node, edge)]
                far.remove(edge)
                path[edge] = init_node
                
    #print(voronoi_cells)
    while considered:
        #print(far)
        #print(considered)
        #print(accepted)
        min_node = None
        for node in considered:
            if min_node is None:
                min_node = node
            elif considered[node] < considered[min_node]:
                min_node = node
        if min_node == None:
            break
        accepted[min_node] = considered[min_node]
        del considered[min_node]
        current_weight = accepted[min_node]
        
        for edge in graph.edges[min_node]:
            try:
                weight = current_weight + graph.distances[(min_node, edge)]
            except:
                weight = current_weight + math.inf
            if (edge not in considered.keys() or weight < considered[edge]) and edge not in accepted.keys() :
                considered[edge] = weight
                path[edge] = min_node
                    
                if edge not in considered:
                    far.remove(edge)
            
        
        
    distances = []
    for i in range(len(accepted)):
        distances.append(accepted[i])
    
    return np.array(distances)