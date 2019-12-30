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
       weight =({'weight' : 10* distance.euclidean(graph[0][edge[0]], graph[0][edge[1]])})
       weighted_graph.append((edge[0],edge[1], weight))
    
    return weighted_graph

def add_connectivity(nodes, weighted_graph, connect):
    nodes_to_connect = np.where(np.isin(nodes[:,0:2], connect[:,0:2])[:,0])[0]
    weight = ({'weight' : 0.01})
    row = 0
    num_rows = np.shape(connect)[0]
    while row < num_rows:
        weighted_graph.append((nodes_to_connect[row], nodes_to_connect[row+1], weight))
        row = row + 2
        
    return weighted_graph

