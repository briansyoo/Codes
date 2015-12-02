# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 13:31:22 2015

A simple 1D brute force cluster identification algorithm. 
This code is intended for learning purposes and has not been optimized. 
Efficiency of this code is of the order O(D*N^2) (D = dimension, N = #points)

@author: byoo1
"""

import numpy as np

#create set of points; user specified
points = [1,2,3,4,4,4,7,8,9,12,13,14]
#distance cutoff for defining a cluster
dcut = 1

#initialize
clusters = []


for i in range(len(points)-1):
    for j in range(i+1,len(points)):

        d= abs(points[i]-points[j])

        
        if d <= dcut:
            #if clusters array is empty (identify first cluster)
            if len(clusters) ==0:            
                if (i and j not in clusters)==0:
                    clusters.append([i,j])
            #if neither i or j are within an already existing cluster        
            elif (i and j not in clusters[-1]):
                if (i not in clusters[-1]):
                    if j not in clusters[-1]:
                        clusters.append([i,j])
                        
            #merge to existing cluster if i or j is already tagged           
            for k in range(len(clusters)):
                
                if (j in clusters[k]) and (i not in clusters[k]):
                    clusters[k].append(i)
                    

                elif (i in clusters[k]) and (j not in clusters[k]):
                    clusters[k].append(j)
                    
                elif (i and j in clusters[k]):
                    continue


#for output purposes
clusters = np.array(clusters)
points =  np.array(points)


if len(clusters)==0:
        print "no clusters identified"
else:
    for i in range(len(clusters)):
        print 'cluster',i,':', points[clusters[i]]


                    

#        