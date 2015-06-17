#!/usr/bin/env python

import ndexClient as nc
import networkx as nx
import re, math, os, sys, operator, random
from scipy.sparse import coo_matrix
import imp
#import kernel_scipy
from optparse import OptionParser
from numpy import genfromtxt, dot
import sys
import math
from array import array
from scipy.sparse.linalg import expm


import kernel_scipy as kernel
import automateHeatKernel as ahk

class WeightedSciPYKernel(kernel.SciPYKernel):
    def __init__(self, network_file, time_T=0.1):
        """ 
        Input:

            network_file - a tab-delimited file in .sif network format:
            <source> <interaction> <target>

        Returns:

            Kernel object.                 

        """

        self.time_T=time_T

        self.labels = {}
        # The number of rows and columns for each kernel
        self.ncols = {}
        self.nrows = {}

        # parse the network, build indexes
        edges, nodes, node_out_degrees,weights = self.parseWeightedNet(network_file)
        num_nodes = len(nodes)
        node_order = list(nodes)
        index2node = {}
        node2index = {}
        for i in range(0, num_nodes):
            index2node[i] = node_order[i]    
            node2index[node_order[i]] = i

        # construct the diagonals
        # SCIPY uses row and column indexes to build the matrix
        # row and columns are just indexes: the data column stores 
        # the actual entries of the matrix
        row = array('i')
        col = array('i')
        data = array('f')
        # build the diagonals, including the out-degree 
        for i in range(0, num_nodes):
            # diag entries: out degree
            degree = 0 
            if index2node[i] in node_out_degrees:
                degree = node_out_degrees[index2node[i]]    
            # append to the end
            # array object: first argument is the index, the second is the data value
            # append the out-degree to the data array
            data.insert(len(data), degree)    
            # build the diagonals
            row.insert(len(row), i)    
            col.insert(len(col), i)

        # add off-diagonal edges 
        for i in range(0, num_nodes):
            for j in range(0, num_nodes):
                if i == j:
                    continue
                if (index2node[i], index2node[j]) not in edges:
                    continue
                # append index to i-th row, j-th column
                row.insert(len(row), i)
                col.insert(len(col), j)
                # -1 for laplacian: i.e. the negative of the adjacency matrix 
                data.insert(len(data), -weights[(index2node[i], index2node[j])])

        # Build the graph laplacian: the CSC matrix provides a sparse matrix format
        # that can be exponentiated efficiently
        L = coo_matrix((data,(row, col)), shape=(num_nodes,num_nodes)).tocsc()
        self.laplacian = L
        self.index2node = index2node
        # this is the matrix exponentiation calculation. 
        # Uses the Pade approximiation for accurate approximation. Computationally expensive.
        # O(n^2), n= # of features, in memory as well. 
        self.kernel = expm(-self.time_T*L)
        self.labels = node_order

    def parseWeightedNet(self, network):
        """
        Parse .sif network, using just the first and third columns
        to build an undirected graph. Store the node out-degrees
        in an index while we're at it. 
        """
        edges = set()
        nodes = set()    
        degrees = {}
        weights = {}
        for line in open(network, 'r'):

            parts = line.rstrip().split("\t")
            source = parts[0]
            weight = float(parts[1])
            target = parts[2]

            # if inputing a multi-graph, skip this
            if (source, target) in edges:
                continue
            if source==target:
                continue

            edges.add((source, target))
            edges.add((target, source))
            nodes.add(source)
            nodes.add(target)

            if source not in degrees:
                degrees[source] = 0
            if target not in degrees:
                degrees[target] = 0

            degrees[source] += weight
            degrees[target] += weight
            weights[(source,target)]=weight
            weights[(target,source)]=weight

        return (edges, nodes, degrees,weights)


def log2plus1(x):
    return(math.log((float(x)+1),2))

def MultiSifToWeightedSif(filename, transform=None, outfile='weighted.sif'):
    f=open(filename,'r')

    edges = {}

    for line in f:
        parts = line.rstrip().split('\t')
        source = parts[0]
        target = parts[2]

        if source ==target:
            continue
        elif ((source, target) or (target,source)) in edges.keys():
            edges[(source,target)]=edges[(source,target)]+1
            edges[(target,source)]=edges[(target,source)]+1
        else:
            edges[(source,target)]=1
            edges[(target,source)]=1

    g=open(outfile, 'w')

    for edge,value in edges.iteritems():
        if not transform:
            g.write(edge[0]+'\t'+str(value)+'\t'+edge[1]+'\n')
        else:
            g.write(edge[0]+'\t'+str(transform(value))+'\t'+edge[1]+'\n')

#main
if __name__=="__main__":
    gene_file='example/some_genes.txt'

    f=open(gene_file,'r')

    get_these=[]

    for line in f:
        get_these.append(line.rstrip())

    ker = WeightedSciPYKernel('example/weighted.sif', time_T=0.1)
    ker.writeKernel('example/weighted_heat_kernel.txt')

    queryVec=ahk.queryVector(get_these,ker.labels)
    diffused=ker.diffuse(queryVec)

    import operator

    sorted_diffused = sorted(diffused.items(), key=operator.itemgetter(1), reverse=True)

    import csv

    writer = csv.writer(open('example/diffused_pancreas_weighted.txt', 'wb'),lineterminator="\n", delimiter='\t')
    for key, value in sorted_diffused:
        writer.writerow([key, value])
