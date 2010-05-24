#!/usr/bin/env python

# Copyright (c) 2010 Ameer Ayoub <ameer.ayoub@gmail.com>

# Attempt to detect graph isomorphy/automorphy through construction of a
# graph's canonical form using node edge degree and BFS tree traversals.

#
# Imports
#

# TODO Make the * imports specific modules/functions
from __future__ import division
from graph_helpers import *
from numpy import *
from time import time
from copy import deepcopy
from copy import copy
from subprocess import Popen
import random
import sys

#
# Variables
#

debug = False                                       # Whether or not to print out debugging messages


#
# Functions
# Refer to graph_helpers module for more function definitions
#

def calc_canonical_form(adjMatrix):
    """Calculates the score lists based on BFS and vertex valence."""
    
    # Calculate the degree of each vertex for reference
    degreeMap = {}
    for i in range(0, adjMatrix.shape[0]):
        sum = 0
        for j in range(0, adjMatrix.shape[1]):
            sum = sum + adjMatrix[i][j]
        degreeMap[i] = sum
    
    # Format of scoreMap is {..., node :[BFS tree], ...}
    # search trees are all initialized to be empty. BFS tree format
    # is a list of lists, where the list index is the level with
    # respect to the tree on which the node was discovered.
    scoreMap = {}
    for key in degreeMap:
        scoreMap[key] = []
    
    if debug:
        print "Initial Map: "
        for n in scoreMap:
            print n, " : ", scoreMap[n]

    
    # Iterate through all the nodes continuing until each node unique BFS tree
    # or a BFS traversal of the tree is complete.
    for v in scoreMap:                              # For each vertex
        
        visited = [v]                                # Initialize the visited node list
        
        fringe = [v]                                 # Initialize the fringe node list
        
        #
        # Start BFS Loop
        #        
        while fringe:
                    
            fringeCopy = copy(fringe)                # Copy the fringe nodes
                    
            fringe = []                                 # Clear the old fringe so we can push new
                                                    # nodes to it
                                                        
            buffer = []
                    
            if debug:                               # Debug
                print "Fringe: ", fringeCopy
                        
            for searchNode in fringeCopy:              
                buffer.append(degreeMap[searchNode])
                for conn in range(0, adjMatrix.shape[1]):
                    # Push adjacent non visited nodes onto fringe
                    if adjMatrix[searchNode][conn] == 1 and not conn in visited:
                        fringe.append(conn)
                        visited.append(conn)     # Mark as visited
            
            buffer.sort()
            scoreMap[v].append(copy(buffer))                  # Update after each iteration
                        
        #        
        # FINISH BFS LOOP
        #
            
    if debug:                                               # Debug
        print "Final Canonical Form: "
        for n in scoreMap:
            print n, ":", scoreMap[n]
    if debug:                                               # Debug seperator
        print "==="
    
    return deepcopy(scoreMap)                               # Copy and return score map

    
def generate_map(adjMatrix, canonicalForm):
    """Calculates all the possible automorphic mappings based on the adjacency
    matrix and the graph's canonical form as calculated by calc_canonical_form."""
    pass

def compare_score_sets(g1, g2):
    """Compares two input canonical form maps checking that they contain the same
    BFS sublists."""
    
    # Remove labels
    g1l = []
    g2l = []
    
    for l1 in g1:
        g1l.append(g1[l1])
    for l2 in g2:
        g2l.append(g2[l2])
        
    cp = deepcopy(g1l)
    
    # Compare
    for li in cp:
        if li in g2l:
            g1l.remove(li)
            g2l.remove(li)
    
    if not g1l and not g2l:
        return True
    if debug:
        print "REMAINING: ", g1l, " \n ", g2l
    return False
        
def generate_map(g, am):
    x = deepcopy(g)
    map = {}
    count = 0
    while x:
        same = []
        min = x.keys()[0]
        # Find Max
        for i in x:
            if x[i] < x[min]:
                min = i
        
        # Find all with the same value
        for i in x:
            if x[i] == x[min]:
                same.append(i)
        
        # Score to order the same values
        samevals = {}
        for s in same:
            # Rather do BFS to find the first values of the nearest connected
            # nodes
            samevals[s] = []
            fringe = [s]
            visited = [s]
            while fringe:
                buffer = []
                fringeCopy = copy(fringe)
                fringe = []
                for searchNode in fringeCopy:
                    if map.has_key(searchNode):
                        buffer.append(map[searchNode])
                    else:
                        buffer.append(-1)
                    for conn in range(0, am.shape[1]):
                        # Push adjacent non visited nodes onto fringe
                        if am[searchNode][conn] == 1 and not conn in visited:
                            fringe.append(conn)
                            visited.append(conn)
                buffer.sort()
                samevals[s].append(copy(buffer))
        
        # Now sort and map
        while samevals:
            # Now sort by samevals and insert into map by that order
            min = samevals.keys()[0]
            for s in samevals:
                if samevals[s] < samevals[min]:
                    min = s
            map[min] = count
            count = count + 1
            del(x[min])
            del(samevals[min])
    return deepcopy(map)

# Test/Debug
if __name__ == "__main__":
    start = time()                                          # Keep track of execution time
    execiso = True                                         # Execute test on isomorphic graphs?
    execrand = False                                         # Execute test on random graphs?

    if execiso:
        max_mult = 1                                            # Real max is 5 when using the graph db
        max_num = 99                                            # Real max is 99 when using the graph db
    
        failed = []                                             # List of tuples of failed graph sets (adj matrices)
        failedId = []                                           # List of the failed graphs to identify the proper one
        results = []                                            # List of boolean results of comparison
        resultsMap = []
        for n in range(1, max_mult+1):
            for i in range(0, max_num+1):
                am1 = read_into_matrix("./graphs/iso_r001_s"+str(20*n)+".A"+str(i).zfill(2))
                m1 = calc_canonical_form(am1)
    
                am2 = read_into_matrix("./graphs/iso_r001_s"+str(20*n)+".B"+str(i).zfill(2))
                m2 = calc_canonical_form(am2)
    
                resultBuff = compare_score_sets(m1, m2)
                results.append(resultBuff)
                
                rm = compare_matrices(am1, generate_map(m1, am1), am2, generate_map(m2, am2))
                resultsMap.append(rm)
                if not rm:
                    print "MAPPING FAILED : ", n*20, " vertex graph : #", i
                    f1 = open("FAILEDA.dot", "w")
                    f1.write(write_to_dot(am1))
                    f1.close()
                    f2 = open("FAILEDB.dot", "w")
                    f2.write(write_to_dot(am2))
                    f2.close()                
                    log_graph("FAILED_"+str(n*20)+"_"+str(i), am1, am2)
                
                if resultBuff == False:
                    failed.append((m1, m2))
                    failedId.append(i)
                    print "NOT ISOMORPHIC : ", i, "."
                
        print "Topology Comparison ("+str(len(results))+" graph sets): "
        print results.count(True)*100/len(results), "% isomorphic."
        print "Mapping Comparison ("+str(len(resultsMap))+" graph sets): "
        print resultsMap.count(True)*100/len(resultsMap), "% mappings identified."
    
    if execrand:
        # Settings
        count = 0
        nodes_min = 10
        nodes_max = 50
        iterations = 1
        
        while count < iterations:
            random.seed()
            nodes = random.randint(nodes_min, nodes_max)
            graph1 = generate_random_graph(nodes)
            graph2 = generate_random_graph(nodes)
            
            score1 = calc_canonical_form(graph1)
            score2 = calc_canonical_form(graph2)
            
            compare = compare_score_sets(score1, score2)
            
            if compare:
                print "ISOMORPHIC GRAPH DETECTED @", count
                log_graph(str(count), graph1, graph2)
                
            count = count + 1
            
    print "Finished in "+str(time()-start)+" sec."