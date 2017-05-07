#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.special


from graphviz import Digraph

curdir = '..'

def visualize_network(A, fname="output", fdir="graph-outputs"):
    dot = Digraph()
    N = A.shape[0]
    for i in range(N):
        dot.node(str(i))
    for i in range(N):
        for j in range(N):
            if A[i,j]:
                dot.edge(str(i), str(j))
    dot.render("{0}/{1}/{2}".format(curdir, fdir, fname))


def load_data():
    curdir = '..'
    with open('{0}/data/training.txt'.format(curdir)) as f:
        #
        data_raw = [ l.split() for l in f]
    keys = data_raw[0]
    data = [ [int(d) for d in l] for l in data_raw[1:] ]
    N = len(data)
    n_vars = len(keys)
    #print(N,n_vars)
    #print(data[0:5])
    data_np = np.array(data)
    if not data_np.shape == (N, n_vars):
        raise ValueError('data in wrong shape')
    return keys, data_np


def initial_network(N=26):
    A_init = np.random.rand(N,N) > 0.5
    # A_init needs to be a DAG
    # first, make it a undirected, no self-loops
    for i in range(N):
        for j in range(N):
            if i > j:
                continue
            if i == j:
                A_init[i,j] = 0
            else:
                A_init[i,j] = A_init[j,i]
    # we have an arbitrary undirected network as a matrix
    # enforce that each node is connected to some node
    for i in range(N):
        if np.all(~A_init[i,:]):
            j = np.random.randint(N)
            A_init[i,j] = 1
            A_init[j,i] = 1


    # make a DAG out of it
    # pick random root
    root = np.random.randint(N)
    list = [root]
    ok = []
    while len(list) > 0:
        j = list.pop(0)
        #print(j, list, ok)
        for i in range(N):
            if i in ok:
                continue
            A_init[i,j] = 0
            list.append(i)
        ok.append(j)

    # TODO I'm quite certain the algorithm produces a correct DAG, but maybe check?

    return A_init


def initial_network_indg(max_indg, N=26):
    """
    graph with in-degree constraint
    """
    pass

    A_init = initial_network(N)

    for i in range(N):
        while np.sum(A_init[:,i]) > max_indg:
            k = np.random.randint(max_indg)
            c = 0
            for t in range(N):
                if A_init[t,i]:
                    c += 1
                if c == (k + 1):
                    A_init[t,i] = 0
                    break

    return A_init


def initial_network_in_outdg(max_indg, max_outdg, N=26):
    """
    graph with in-degree constraint
    """
    pass

    #todo
    A_init = initial_network(N)

    for i in range(N):
        while np.sum(A_init[:,i]) > max_indg:
            k = np.random.randint(max_indg)
            c = 0
            for t in range(N):
                if A_init[t,i]:
                    c += 1
                if c == (k + 1):
                    A_init[t,i] = 0
                    break

    for i in range(N):
        while np.sum(A_init[i,:]) > max_outdg:
            k = np.random.randint(max_outdg)
            c = 0
            for t in range(N):
                if A_init[i,t]:
                    c += 1
                if c == (k + 1):
                    A_init[i,t] = 0
                    break

    return A_init


def check_dagness(A, start, added):
    """
    Return true if A still a DAG after adding i->j. Start search from (the recently added) connection i->j
    """
    list = [added]
    N = A.shape[0]
    while len(list) > 0:
        current = list.pop(0)
        for i in range(N):
            if A[current, i]:
                if i == start:
                    return False
                list.append(i)

    return True
