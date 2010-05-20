#!/usr/bin/env python

# Copyright (c) 2010 Ameer Ayoub <ameer.ayoub@gmail.com>

# Attempt to detect graph isomorphy/automorphy through construction of a
# graph's canonical form using node edge degree and BFS tree traversals.

#
# Imports
#
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
#
def calc_score_map(adjMatrix):
    """Calculates the standard label map based on node weights assigned by
    node degree."""
    
    # Calculate the degree of each node, this is the degree class
    weightMap = {}
    for i in range(0, adjMatrix.shape[0]):
        sum = 0
        for j in range(0, adjMatrix.shape[1]):
            sum = sum + adjMatrix[i][j]
        weightMap[i] = sum
    
    # Format of weightMapRev is (node, score)
    # scores are all initialized to zero
    weightMapRev = {}
    for key in weightMap:
        if weightMapRev.has_key(weightMap[key]):
            weightMapRev[weightMap[key]][key] = 0
        else:
            weightMapRev[weightMap[key]] = {key: 0}
    if debug:
        print "Initial Map: "
        for n in weightMapRev:
            print n, " : ", weightMapRev[n]
    
    maxScore = 0
    for n in weightMap:
        maxScore = maxScore + weightMap[n]
    if debug:
        print "MAX: ", maxScore
    
    # Iterate through all the nodes of the same weight, continuing until
    # each node in the same weight class has a unique weight value or a BFS
    # traversal of the tree is complete.
    for k in weightMapRev:                          # For each weight class
        
        visited = {}                                # Initialize the visited node list
                                                    # for each node
        for n in weightMapRev[k]:
            visited[n] = [n]
        
        fringe = {}                                 # Initialize the fringe node list
                                                    # for each node
        for n in weightMapRev[k]:
            fringe[n] = [n]
        
        roundCounter = 0                            # Round counter to score nodes
                                                    # that finish in a group together
        
        # Continue execution until we have traveresed all trees or until
        # all the nodes have a unique score.
        #
        # STARTS BFS LOOP
        #
        while not allUnique(weightMapRev[k]) and not nonUniqueSublistsAreEmpty(fringe, weightMapRev[k]):
        
            # Selection process for picking the nodes that are to execute a round
            # of BFS. BFS is executed only on nodes of duplicate score, scores are
            # picked with preference to the largest score.
            
            processList = []                        # List of nodes to run BFS on
            
            # Pick initial value
            for m in weightMapRev[k]:
                if not isUnique(m, weightMapRev[k]) and fringe[m]:
                    pick = m

            # Pick largest one
            for m in weightMapRev[k]:
                if weightMapRev[k][m] >= weightMapRev[k][pick] and not isUnique(m, weightMapRev[k]) and fringe[m]:
                    pick = m
                
            if debug:                                       # Debug
                print "Picked Value: ", weightMapRev[k][pick]
            
            # Add all the nodes with the picked score onto the list to be processed
            for m in weightMapRev[k]:
                if weightMapRev[k][m] == weightMapRev[k][pick]:
                    processList.append(m)
                    
            if debug:                                       # Debug
                print "Process List: ", processList
                    
            for m in processList:                           # For all we should calculate BFS for
                if fringe[m]:                               # If the node rooted BFS traversal isn't complete
                    
                    if debug:                               # Debug
                        print "BFS On Node ", m
                    
                    sum = weightMapRev[k][m]                # Reset this nodes sum to its previous score
                    
                    fringeCopy = copy(fringe[m])            # Copy the fringe nodes
                    
                    fringe[m] = []                          # Clear the old fringe so we can push new
                                                            # nodes to it
                    
                    if debug:                               # Debug
                        print "Fringe: ", fringeCopy
                        
                    for searchNode in fringeCopy:           # Do BFS
                        
                        if debug:                           # Debug
                            print weightMap[searchNode], " + ",
                            
                        sum = sum + weightMap[searchNode]   # Count towards sum
                        
                        for conn in range(0, adjMatrix.shape[1]):
                            # Push adjacent non visited nodes onto fringe
                            if adjMatrix[searchNode][conn] == 1 and not conn in visited[m]:
                                fringe[m].append(conn)
                                visited[m].append(conn)     # Mark as visited
                                
                    weightMapRev[k][m] = sum                # Update weight map
                    if debug:                               # Debug
                        print " = ",sum
            if debug:                                       # Debug
                print weightMapRev[k]
            
            # Adjust score of nodes that have finished to ensure that they are done
            for m in weightMapRev[k]:
                if weightMapRev[k][m] == maxScore:
                    weightMapRev[k][m] = maxScore + roundCounter
                    
            roundCounter = roundCounter + 1                 # Update the round counter
            
            #        
            # FINISH BFS LOOP
            #
            
        # Relabel symmetric groups by group cardinality
        groups = []                                         # List to hold groups
        
        for l1 in weightMapRev[k]:
            if not isUnique(l1, weightMapRev[k]) and weightMapRev[k][l1] not in groups:
                groups.append(weightMapRev[k][l1])          # Add non-unique values to groups list
                
        if debug:
            print "Groups: ", groups                        # Debug
            
        # Sanity Check: All nodes in a group connect to the same node(s)
        # Note this only would work if all the symmetric groups are connected
        # this doesn't work in the 'nested' symmetric subtree.
        """
        flag = True
        for g in groups:
            group_fringe_buffer = []
            for n in weightMapRev[k]:
                if weightMapRev[k][n] == g:
                    if not group_fringe_buffer:
                        # Push adjacent nodes into the buffer
                        for conn in range(0, adjMatrix.shape[0]):
                            if adjMatrix[n][conn] == 1:
                                group_fringe_buffer.append(conn)
                    else:
                        # Check to make sure nodes are the same
                        for conn in range(0, adjMatrix.shape[0]):
                            if adjMatrix[n][conn] == 1:
                                if conn not in group_fringe_buffer:
                                    print "FAIL: Group fringe does not contains node ", conn
                                    flag = False
                            if adjMatrix[n][conn] == 0:
                                if conn in group_fringe_buffer:
                                    print "FAIL: Group fringe contains node ", conn
                                    flag = False
        if flag and debug:
            print "PASS: All groups are symmetric proceeding with group scoring."
        """                            
        
        # Groups now contains all the symmetric node group scores, now lets relabel them
        # Relabeling occurs through similar iterative BFS rooted at the group
        
        """
        weightMapBufferGroup = {}
        for g in groups:
            for n in weightMapRev[k]:
                weightMapBufferGroup[n] = 0
        
        fringeGroup = {}
        visitedGroup = {}
        for g in groups:
            fringeGroup[g]
        for n in weightMapRev[k]:
            if weightMapRev[k][n] == g:
                fringeGroup[g].append(n)
                visitedGrou[g].append(n)
        """
        
        # OLD
        for g in groups:                                    # For each group
            
            weightMapCop = weightMapRev[k].copy()           # Initalize a buffer of the score map
                                                            # to copy the new scores to
                                                            
            count = 0                                       # Count the instances of the score
            for l1 in weightMapRev[k]:
                if weightMapRev[k][l1] == g:
                    count = count + 1
            
            for l1 in weightMapRev[k]:                      # Relabel similar scored nodes
                if weightMapRev[k][l1] == g:
                    weightMapCop[l1] = maxScore + count     # Relabel to the buffer to avoid
                                                            # relabelling multiple times
                    
            weightMapRev[k] = weightMapCop.copy()           # Copy buffer to main score map
            
        if debug:                                           # Debug seperator
            print "---"            
    
    if debug:                                               # Debug
        print "Final Weight Map: "
        for n in weightMapRev:
            print n, ":", weightMapRev[n]
    if debug:                                               # Debug seperator
        print "==="
    
    returnMap = deepcopy(weightMapRev)                      # Copy the return score map
    
    nodeMap = {}
    counter = 0
    for i in sorted(weightMapRev):
        while weightMapRev[i]:
            max = weightMapRev[i].keys()[0]
            for j in weightMapRev[i]:
                if weightMapRev[i][j] > weightMapRev[i][max]:
                    max = j
            nodeMap[max] = counter
            del weightMapRev[i][max]
            counter = counter + 1
            
    return (returnMap, nodeMap)

def compare_score_sets(weightMap1, weightMap2):
    """Compares two input weight maps checking that the score structure
    is the same."""
    
    weightMap2Copy = deepcopy(weightMap2)
    for k in weightMap1:
        buffer1 = []
        buffer2 = []
        if not weightMap2.has_key(k):                       # Buffer 1 has a node set buffer 2 doesn't
            if debug:
                print "Weight map 1 contains unmatched node set(s)."
            return False
        for l1 in weightMap1[k]:
            buffer1.append(weightMap1[k][l1])               # Copy all elements to buffer1 from weightMap1
        for l2 in weightMap2[k]:
            buffer2.append(weightMap2[k][l2])               # Copy all elements to buffer2 from weightMap2
        if debug:
            print buffer1, " X ", buffer2
        for s in buffer1:                                   # Remove corresponding elements one by one from both buffers
            if s in buffer2:
                buffer2.remove(s)
            else:                                           # Buffer 2 doesn't have an element that weightMap 1 did
                if debug:
                    print "Weight map 1 still contains unmatched element(s)."
                    print buffer1, " X ", buffer2
                return False
        if buffer2:                                         # Buffer 2 contains more elements that buffer 1 didn't have
            if debug:                                       # Debug
                print "Weight map 2 still contains unmatched element(s)."
                print buffer1, " X ", buffer2
            return False
        del(weightMap2Copy[k])                               # If used all elements then delete
    if not weightMap2Copy:                                   # If used all node sets in the weight map then we completed a match
        return True
    else:
        if debug:
            print "Weight map 2 still contains unmatched node set(s)."
            print weightMap2Copy
        return False


# Test/Debug
if __name__ == "__main__":
    start = time()                                          # Keep track of execution time

    """
    max_mult = 1                                            # Real max is 5 when using the graph db
    max_num = 0                                             # Real max is 99 when using the graph db
    
    failed = []                                             # List of tuples of failed graph sets (adj matrices)
    failedId = []                                           # List of the failed graphs to identify the proper one
    results = []                                            # List of boolean results of comparison
    resultsMap = []
    for n in range(1, max_mult+1):
        for i in range(41, 42):
            am1 = read_into_matrix("C:/Users/Ameer/My Documents/Komodo Edit Projects/GraphRelabeling/graphs/iso_r001_s"+str(20*n)+".A"+str(i).zfill(2))
            m1 = calc_score_map(am1)
            wm1 = m1[0]
            nm1 = m1[1]

            am2 = read_into_matrix("C:/Users/Ameer/My Documents/Komodo Edit Projects/GraphRelabeling/graphs/iso_r001_s"+str(20*n)+".B"+str(i).zfill(2))
            m2 = calc_score_map(am2)
            wm2 = m2[0]
            nm2 = m2[1]

            resultBuff = compare_score_sets(wm1, wm2)
            results.append(resultBuff)
            resultsMapBuff = compare_matrices(am1, nm1, am2, nm2)
            resultsMap.append(resultsMapBuff)
            
            if resultBuff == False:
                failed.append((am1, am2))
                failedId.append(i)
                count1 = 0
                for q in range(0, am1.shape[0]):
                    for k in range(0, am1.shape[1]):
                        if am1[q][k] == 1:
                            count1 = count1 +1
                count2 = 0
                for q in range(0, am2.shape[0]):
                    for k in range(0, am2.shape[1]):
                        if am2[q][k] == 1:
                            count2 = count2 +1
                print "NOT ISOMORPHIC : ", i, "; ", count1, "edges vs ", count2, "edges."
            
    print "Small Topology Comparison ("+str(len(results))+" graph sets): "
    print results.count(True)*100/len(results), "% isomorphic."
    print "Small Mapping Comparison ("+str(len(results))+" graph sets): "
    print resultsMap.count(True)*100/len(resultsMap), "% mappings identified."
    
    if failed:
        minimap = {0: "A", 1: "B"}
        print "Writing failed graphs to file..."
        count = 0
        for s in range(0, len(failed)):
            for g in minimap:
                filepath = "./failed_out/"+str(failedId[count]).zfill(2)+str(minimap[g])+".dot"
                print "Writing %s" % filepath
                f = open(filepath, 'w')
                f.write(write_to_dot(failed[s][g]))
                f.close()
            count = count + 1
    """
    
    # Settings
    count = 0
    nodes_min = 10
    nodes_max = 50
    iterations = 1
    
    while count < iterations:
        random.seed()
        nodes = random.randint(nodes_min, nodes_max)
        graph1 = generate_random_graph(nodes)
        #print graph1
        graph2 = generate_random_graph(nodes)
        #print graph2
        score1 = calc_score_map(graph1)[0]
        score2 = calc_score_map(graph2)[0]
        compare = compare_score_sets(score1, score2)
        if compare:
            print "ISOMORPHIC GRAPH DETECTED @", count
            filepath1 = "./iso_graphs/"+str(count)+"A.matrix"
            f1 = open(filepath1, 'w')
            f1.write(graph1)
            f1.close
            filepath2 = "./iso_graphs/"+str(count)+"B.matrix"
            f2 = open(filepath2, 'w')
            f2.write(graph2)
            f2.close
            filepath3 = "./iso_graphs/"+str(count)+"A.graphdb"
            f3 = open(filepath3, 'w')
            f3.write(write_to_vflib(graph1))
            f3.close
            filepath4 = "./iso_graphs/"+str(count)+"B.graphdb"
            f4 = open(filepath4, 'w')
            f4.write(write_to_vflib(graph2))
            f4.close
        count = count + 1
        if count%10 == 0:
            print count/iterations, "% complete"
            sys.stdout.flush()
    
    print "Finished in "+str(time()-start)+" sec."