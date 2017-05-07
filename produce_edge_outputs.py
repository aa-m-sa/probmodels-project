#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

curdir = '..'

abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

USE_SCALES = True

if __name__ == "__main__":
    count_connections = dict()
    count_connections_ind = dict()
    #for outletter in "ABCDEFGH":

    sworst = None
    sbest = None
    for outletter in "ABCDEFGHIJKLMNOPQRST":
        for series in ['bdeu']:
            outname = "{0}_newfix_2500_{1}".format(series, outletter)
            #f = np.load(curdir+'/testoutput-2017-04-12b/testrun{0}.npz'.format(outname))
            f = np.load(curdir+'/testoutput-2017-04-26/testrun_{0}.npz'.format(outname))

            s = np.max(f['s_history'])
            if not sworst or sworst > s:
                sworst = s

            if not sbest or sbest < s:
                sbest = s

            f.close()


    for outletter in "ABCDEFGHIJKLMNOPQRST":
        #for series in ['bic', 'aic']:
        for series in ['bdeu']:
            #outname = "{0}_dnp_2500_{1}".format(series, outletter)
            outname = "{0}_newfix_2500_{1}".format(series, outletter)
            #f = np.load(curdir+'/testoutput-2017-04-12b/testrun{0}.npz'.format(outname))
            f = np.load(curdir+'/testoutput-2017-04-26/testrun_{0}.npz'.format(outname))

            A_best = f['A_best']
            s = np.max(f['s_history'])
            scale = np.abs(sbest) / np.abs(s)

            for i, a in enumerate(abc):
                for j, b in enumerate(abc):
                    if i == j:
                        continue
                    if (a,b) not in count_connections:
                        count_connections[(a,b)] = 0
                        count_connections_ind[(i,j)] = 0

                    k = 1
                    if USE_SCALES:
                        k = k * scale
                    if A_best[i,j]:
                        count_connections[(a,b)] += k
                        count_connections_ind[(i,j)] += k

    sorted_counts = sorted(count_connections.items(), key= lambda x: x[1], reverse=True)
    sorted_counts_ind = sorted(count_connections_ind.items(), key= lambda x: x[1], reverse=True)
    At = np.zeros((len(abc), len(abc)))
    for edge, c in sorted_counts:
        #print(edge[0], edge[1], c)
        print(edge[0], edge[1])
    for edge, c in sorted_counts_ind:
        if c > 1:
            At[edge[0],edge[1]] = 1

    #np.savez(curdir+'/testoutput-2017-04-12b/thres_ge2_aicbic_average_16_2500.npz', At = At)
    #np.savez(curdir+'/testoutput-2017-04-12b/thres_ge2_aic_average_8_2500.npz', At = At)
