#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

curdir = '..'

abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


if __name__ == "__main__":
    fname = 'out-long-mi7.txt'

    conn_probs = dict()

    with open('{0}/disco-output/{1}'.format(curdir, fname)) as f:
        #
        data_raw = [ l.split() for l in f]

    for l in data_raw:
        node_from = int(l[0])
        node_to = int(l[2])

        value = np.double(l[3])

        conn_probs[(node_from, node_to)]  = value

    sorted_probs = sorted(conn_probs.items(), key = lambda x: x[1], reverse=True)

    A = np.zeros((len(abc),len(abc)))

    for k in sorted_probs:
        node_from = k[0][0]
        node_to = k[0][1]
        print(abc[node_from], abc[node_to])
        #print(abc[node_from], abc[node_to], k[1])
        value = k[1]
        if value > 0.4:
            if A[node_to, node_from] == 1.0:
                continue
            A[node_from, node_to] = 1.0


    outname = '{0}_04'.format(fname)
    np.savez('{0}/disco-output/{1}'.format(curdir, outname), A=A)
