# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 13:31:22 2015

A simple 1D brute force single cluster identification algorithm.
The algorithm will identify all neighboring points recursively
for a user-provided target (alive) point. 
 
This code is intended for learning purposes and has not been optimized. 

@author: byoo1
"""

import numpy as np

#create set of points; user specified
points = [1,2,3,4,4,4,4,7,8,9,12,13,14]

dcut = 1
#target point index
alive = int(np.floor(np.random.randint(len(points))))


#initialize
clusters = []

def find_neighbors(alive,points):
    #reset neighbor list
    neighbors=[]
    for i in range(len(points)):
        d = abs(points[i]-points[alive])
        if d <=dcut and i != alive:
            if i not in clusters:
                clusters.append(i)
                neighbors.append(i)
        #once end of loop is reached
        if i ==(len(points)-1):
            #loop through neighbor indices recursively
            for j in neighbors:
                find_neighbors(j,points)        
        else:
            continue
        
    return

find_neighbors(alive,points)
clusters = np.array(clusters)
points = np.array(points)

#for output
print "alive index:", alive
print "alive value:", points[alive]
print "points:"
print "cutoff distance:", dcut
for i in range(len(points)):
    print points[i]


cluster_points = points[clusters]
print "cluster points:"
for i in range(len(cluster_points)):
    print cluster_points[i]





