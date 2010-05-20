#!/usr/bin/env python

# Copyright (c) 2010 Ameer Ayoub <ameer.ayoub@gmail.com>

# Helper functions for working with graphs.

#
# Imports
#
from __future__ import division
from numpy import *
from struct import *
from copy import deepcopy
from copy import copy
from pygraph.classes.graph import graph
from pygraph.readwrite.dot import write
import random

#
# Functions
#

#
# Files and Formats
#

def read_short(file):
    """Reads in a short from a file like object, returns as a python integer."""
    
    return copy(unpack('h', file.read(2))[0])       # Read in two bytes as short
                                                    # See http://docs.python.org/library/struct.html
                                                    

def read_into_matrix(fname):
    """Reads a graph database file and returns a corresponding adjacency
    matrix (of type numerical array). Format of the input file is described
    in detail in the graph database documentation
    http://amalfi.dis.unina.it/graph/db/doc/graphdb-3.html the code here
    is largely based on the c example code found in the graph database
    documentation."""
    
    f = open(fname, 'rb')                           # Open file in binary reading format
    
    nodes = read_short(f)                           # Read in the number of nodes
                                                    
    adjMatrix = zeros((nodes, nodes), dtype=int)    # Initalize a 2D array of size nodes x nodes
                                                    # This will be the adjacency matrix
    
    for i in range(0, nodes):                       # For each node
        edges = read_short(f)                       # Read in the number of edges
        for j in range(0, edges):                   # For each edge
            dest = read_short(f)                    # Read the destination node
            adjMatrix[i][dest] = 1;                 # Mark the connection in the adjacency matrix
            adjMatrix[dest][i] = 1;                 # Mark the symmetric connection for ease of lookup
            
    f.close()                                       # Close the file
    return adjMatrix


def write_to_vflib(adj_matrix):
    """Returns the appropriate VFLib file format from an adjacency matrix."""
    
    vflib_string = pack('h', adj_matrix.shape[0])
    for i in range(0, adj_matrix.shape[0]):
        edge_buffer = []
        count = 0
        for j in range(0, i+1):
            if adj_matrix[i][j] == 1:
                count = count + 1
                edge_buffer.append(j)
        vflib_string = vflib_string + pack('h', count)
        for edge in edge_buffer:
            vflib_string = vflib_string + pack('h', count)
    return vflib_string


def write_to_dot(adj_matrix):
    """Returns the dot file format string from an adjacency matrix for use with
    graphviz."""
    
    gr = graph()
    gr.add_nodes(range(0, adj_matrix.shape[0]))
    for i in range(0, adj_matrix.shape[0]):
        for j in range(0, i+1):
            if adj_matrix[i][j] == 1:
                gr.add_edge([i, j])
    return write(gr)

 
#
# Graph Generators
#

def generate_random_graph(nodes, freq_inv=2):
    """Returns a graph containing nodes number of nodes with randomly selected
    connections, every node is approximately  connected to 50% of all other
    nodes, making the graphs very densely connected by default. To set the
    frequency of connections pass in connections frequency inverse as freq_inv
    e.g. 2 is 1/2 frequency."""
    
    random.seed()                                   # Seed the random generator from time
    return_matrix = zeros((nodes, nodes), dtype=int)# Initialize the adjacency matrix to zeroes
    for i in range(0, return_matrix.shape[0]):      # Loop through every connection
        for j in range(0, i+1):
            if random.randint(0, freq_inv-1)%freq_inv == 0:
                conn = 1                            # If the random number mod freq_inv is 0
                                                    # then set a connection (this occurs with
                                                    # 1/freq_inv frequency.)
                
                return_matrix[i][j] = conn          # Only set connection if its 1 otherwise
                return_matrix[j][i] = conn          # it's already 0. 
    return return_matrix


#
# Comparison Functions
#

def compare_matrices(adjMatrix1, nodeMap1, adjMatrix2, nodeMap2):
    """Takes in two adjacency matrices and two node mapping dictionaries and
    checks for equivalency. The matrices are unlabeled numerical arrays and
    the mapping dictionaries are standard python dictionaries containing
    mappings to use when comparing the matrices, e.g. {..., 1: 2, ...} implies
    that the node in column 1 is actually to be refered to as node 2 when
    comparing the matrices. Returns true if the two matrices are equivalent
    under the specified node mapping, false otherwise."""
    
    if adjMatrix1.shape != adjMatrix2.shape:        # Check that the matrices are the same dimensions
        print "Graphs not the same size."
        return False
    
                                                    # Reverse nodeMap2 for reverse translation
    nodeMap2_rev = {}                               # Initialize nodeMap2_rev to empty dictionary
    for key in nodeMap2:                            # Iterate over all the keys
        nodeMap2_rev[nodeMap2[key]] = key           # Create a new entry in the dictionary {..., value:key}
    
    for i in range(0, adjMatrix1.shape[0]):         # Loop through all the columns
        for j in range(0, adjMatrix1.shape[1]):     # Loop through all the rows
            
            r = nodeMap1[i]                         # Resolve column mapping for matrix1
            c = nodeMap1[j]                         # Resolve row mapping for matrix1
            
            # Reverse lookup the col/row mapping on matrix2 and compare connections
            if adjMatrix1[i][j] != adjMatrix2[nodeMap2_rev[r]][nodeMap2_rev[c]]:
                return False                        # If they aren't equal then return false
            
    return True                                     # Otherwise if we looped through everything
                                                    # and no edges are different then return
                                                    # true


#
# Computational Helpers
#

def isUnique(key, dic):
    """Returns True if the value of the key parameter passed in is unique within
    the dictionary parameter."""
    
    try:
        for d in dic:
            if dic[key] == dic[d] and d != key:
                return False
        return True
    except KeyError:
        print "Invalid dictionary key \"%s\"when checking uniqueness. Check usage" % key


def allUnique(dic):
    """Returns True if all the values in the dictionary are unique."""
    
    for d in dic:
        if not isUnique(d, dic):
            return False
    return True


def sublistsAreEmpty(lis):
    """Checks whether all sublists of a collection are empty. Not recursive."""
    
    for l in lis:
        if lis[l]:
            return False
    return True


def nonUniqueSublistsAreEmpty(lis, dic):
    """Checks whether all the sublists belonging to keys with nonunique values
    are exhausted."""
    
    try:
        for d in dic:
            if not isUnique(d, dic):
                if lis[d]:
                    return False
        return True
    except KeyError:
        print "Sublist check failed, ensure dictionary/list coordination."
