
# coding: utf-8

# In[1]:

import numpy as np
#import dateutil.parser as dparser
#import datetime
import matplotlib.pyplot as plt
#import urllib2
import json
#from collections import Counter
#import pprint
#import time
import os, sys
import csv
import scipy, scipy.stats
#import random
#import re
#import networkx as nx
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import linkage, dendrogram
from itertools import groupby
#matplotlib inline 


# In[8]:

All_100_null_batches =json.load(open('temp_stats/Biased_Random_Walks_All_100.txt'))


# In[2]:

def get_n_clusters_winner_matters(n):
    T = sch.fcluster(link_hierarchy_ward_winner_matters, n, 'maxclust')
    
    # calculate labels
    labels=list('' for i in range(1310))
    for i in range(1310):
        labels[i]=(i,T[i])
    
    # depth_to
    # cluster[d] = labels
    
    clusters = []
    for i in range(n):
        clusters.append([])

    for x in labels:
        clusters[x[1]-1].append(x[0])
        
    return clusters


# In[3]:

def avg_intra_cluster(cluster, dist_matrix):
    d = []
    for x in cluster:
        for y in cluster:
            if x!=y:
                d.append(dist_matrix[x][y])
    return np.mean(d)


# In[4]:

def how_many_clusters(worm_set,index):
    
    #first measure the distances bewtween them all
    
    biased_normed_dist_winner_matters = np.zeros((1310,1310))
    done_list = set()
    for i in range(1310):
        worm1 = worm_set[i]
        if i%600==0:
            print i
        for j in range(1310):
            if (i,j) not in done_list:
                worm2 = worm_set[j]
                d = compare_worms_normed(worm1,worm2,orient = 1)
                biased_normed_dist_winner_matters[i][j] = d
                biased_normed_dist_winner_matters[j][i] = d
                done_list.add((i,j))
                done_list.add((j,i))
            
                
    #print ' done computing distances'
    #save
    np.savetxt("temp_stats/null_game_dists/Null_Game_Distances_"+str(index)+".csv", biased_normed_dist_winner_matters, delimiter=",")
                
        
        
    # convert to array
    #biased_normed_dist_winner_matters_array = biased_normed_dist_winner_matters
    
    biased_normed_dist_winner_matters[0][0] = 0.0
    link_hierarchy_ward_winner_matters = linkage(np.array(biased_normed_dist_winner_matters),method = 'ward')
    
    #now compute the clusters
    
    number_of_motifs = {}
    flag_11 = False
    flag_9 = False
    d= 1

    while not flag_11 or not flag_9:

        T = sch.fcluster(link_hierarchy_ward_winner_matters, d, 'maxclust')

        # calculate labels
        labels=list('' for i in range(1310))
        for i in range(1310):
            labels[i]=(i,T[i])

        #depth_to_cluster[d] = labels

        clusters = []
        for i in range(d):
            clusters.append([])

        for x in labels:
            clusters[x[1]-1].append(x[0])


        i_dists = []
        for cluster in clusters:
            i_dists.append(avg_intra_cluster(cluster,biased_normed_dist_winner_matters ))

        avg_pts = np.nanmean(i_dists)
        if avg_pts<11 and not flag_11:
            number_of_motifs[11]=d
            flag_11 = True


        if avg_pts<9 and not flag_9:
            number_of_motifs[9]=d
            flag_9 = True

        d+=1
    return number_of_motifs



# In[5]:

def compare_worms_normed(w1,w2,orient = None):
    len1 = len(w1)
    len2 = len(w2)
    
    if len1>len2:
        reference = w1
        compare = w2
    else:
        reference = w2
        compare = w1
        
    compare = compare+[compare[-1]]*(abs(len1-len2)) # make worms the same length by appending on to the shorter one
    
    # flip up the reference worm so that it is oriented to the winning team on top
    if reference[-1]<0:
        reference = [-x for x in reference]
    
    # option to orient the compare worm in the same way as the reference
    if orient:
        if compare[-1]<0:
            compare = [-x for x in compare]
    
    diff1 = 0 # difference in scores
    diff2 = 0
    
    for t in range(len(reference)):
        diff1+=abs(reference[t]-compare[t])
        diff2+=abs(reference[t]+compare[t]) # flip up the compare worm
        
    dist= min([diff1,diff2])
    if orient:
        dist = diff1
    
    normed_dist = float(dist)/max([len1,len2]) # normalize by the length of the reference game
    
    return normed_dist


# In[6]:

number_of_motifs_11 = [14,
 10,
 13,
 13,
 11,
 10,
 12,
 15,
 12,
 11,
 12,
 11,
 14,
 15,
 12,
 14,
 11,
 13,
 13,
 12,
 9,
 11,
 9,
 15,
 12,
 15,
 14,
 10,
 13,
 14,
 12,
 13,
 11,
 10]

number_of_motifs_9 = [65,
 45,
 57,
 49,
 45,
 46,
 56,
 56,
 45,
 48,
 52,
 45,
 44,
 57,
 49,
 41,
 54,
 57,
 48,
 63,
 44,
 58,
 41,
 55,
 47,
 59,
 48,
 44,
 50,
 49,
 51,
 55,
 45,
 40]


# In[7]:

# start_at_index = 34
start_at_index = int(sys.argv[1])

# In[ ]:

ct = start_at_index
for worm_set in All_100_null_batches[start_at_index:]:
    motifs_dict = how_many_clusters(worm_set,ct)
    
    number_of_motifs_11.append(motifs_dict[11])
    number_of_motifs_9.append(motifs_dict[9])
    
    #save the progress for goodness sake!
    with open('temp_stats/number_of_motifs_11_pts', 'w') as outfile:
        json.dump(number_of_motifs_11, outfile)
    
    with open('temp_stats/number_of_motifs_9_pts', 'w') as outfile:
        json.dump(number_of_motifs_9, outfile)
    
#     print
#     print motifs_dict
    
#     print ct, ' batches done'
#     print 
    ct+=1


# In[ ]:



