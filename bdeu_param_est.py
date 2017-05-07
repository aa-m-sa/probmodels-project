#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from score_structure_learn import load_data
from datetime import datetime

import warnings

import stuff

curdir = '..'

# baysian MAP parameter estimator with BDEU prior for a given network

def params_estimate_bayesian_bdeu(A, alpha, data):
    #
    n_parents = np.sum(A,0)

    qi_terms = np.power(3,n_parents)
    ri = 3 # constant for all i

    N = A.shape[0] # note, number of variables
    data_len = data.shape[0]

    alpha = np.double(alpha)

    est = dict()
    for i in range(N):
        parents = np.where(A[:,i])[0]
        confs_seen = dict() # N_ij
        confs_seen_by_ival = dict() # N_ijk
        for n in range(data_len):
            ival = data[n,i]

            conf_orig = data[n,:][parents]
            conf = conf_orig.tostring() # 'hash'

            if conf not in confs_seen:
                confs_seen[conf] = 1
            else:
                confs_seen[conf] += 1

            if (conf, ival) not in confs_seen_by_ival:
                confs_seen_by_ival[(conf, ival)] = 1
            else:
                confs_seen_by_ival[(conf, ival)] += 1

        tmp = 0.0
        for conf in confs_seen:
            Nij = confs_seen[conf]
            alpha_ij = alpha/qi_terms[i]
            alpha_ijk = alpha_ij/ri # ri constant for all i
            for ival in range(3):
                if (conf, ival) not in confs_seen_by_ival:
                    Nijk = np.finfo(np.double).eps
                else:
                    Nijk = confs_seen_by_ival[(conf, ival)]
                Nijk = np.double(Nijk)
                est_ijk = (Nijk + alpha_ijk)/(Nij + alpha_ij)
                # see slides for Bayesian predictive
                # todo check why -1 and -ri terms would be missing?
                est[(i,conf,ival)] = est_ijk
                # note that there is no value set for configurations not seen in tr data!
                # should be (alpha_ijk - 1) / (alpha_ij - ri)

    return est


def predictive_bayesian_bdeu(A, alpha, predict, tr_data):
    # predict is a vector of variable configurations, task is to calculate a  predictive distribution on them
    #
    n_parents = np.sum(A,0)

    qi_terms = np.power(3,n_parents)
    ri = 3 # constant for all i

    N = A.shape[0] # note, number of variables

    alpha = np.double(alpha)
    est = params_estimate_bayesian_bdeu(A, alpha, tr_data)
    d_preds = []
    n_veryweird = 0
    n_allweird = 0
    for d in predict:
        log_tmp = 0.0
        for i,v in enumerate(d):
            parents_i = np.where(A[:,i])[0]
            parents_st = d[parents_i.astype(int)] # conf
            conf_id = parents_st.tostring()
            was_in_est = False
            if (i,conf_id,v) in est:
                p = est[(i,conf_id,v)]
                was_in_est = True
            else:
                alpha_ij = alpha/qi_terms[i]
                alpha_ijk = alpha_ij/ri # ri constant for all i
                #p  = (alpha_ijk - 1)/(alpha_ij - ri)
                p = alpha_ijk/alpha_ij

            if p <= 0.0:
                #print('Weird p')
                #print(p)
                #print(d)
                #print(i,v)
                #print(parents_st)
                #print(was_in_est)
                if was_in_est:
                    n_veryweird += 1
                n_allweird +=1

            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                try:
                    log_tmp += np.log(p) # p = p(i,v | parents_i))
                except Warning:
                    #print('RuntimeWarning')
                    #print(p, log_tmp)
                    pass
        d_preds.append(np.exp(log_tmp))

    #print('final weird')
    #print(n_veryweird, n_allweird)
    return d_preds, est


def process_python_ensembles():
    #f = np.load(curdir+'/testoutput/thres_ge2_aicbic_average_12_2300.npz')
    #A = f['At']
    USE_SCALES = True

    keys, tr_data = load_data()

    with open('{0}/data/test.txt'.format(curdir)) as f:
        #
        data_raw = [ l.split() for l in f]
    pred_data = [ [int(d) for d in l] for l in data_raw[1:] ]
    pred_data = np.array(pred_data)

    all_preds = np.zeros(pred_data.shape[0])

    sworst = None
    sbest = None
    for outletter in "ABCDEFGHIJKLMNOPQRST":
        for series in ['bdeu']:
            #print(outletter)
            outname = "{0}_newfix_2500_{1}".format(series, outletter)
            #f = np.load(curdir+'/testoutput-2017-04-12b/testrun{0}.npz'.format(outname))
            f = np.load(curdir+'/testoutput-2017-04-26/testrun_{0}.npz'.format(outname))

            s = np.max(f['s_history'])
            if not sworst or sworst > s:
                sworst = s

            if not sbest or sbest < s:
                sbest = s

            f.close()

    sum_scales = 0

    for outletter in "ABCDEFGHIJKLMNOPQRST":
        for series in ['bdeu']:
            #print(outletter)
            #print(str(datetime.now()))
            outname = "{0}_newfix_2500_{1}".format(series, outletter)
            #f = np.load(curdir+'/testoutput-2017-04-12b/testrun{0}.npz'.format(outname))
            f = np.load(curdir+'/testoutput-2017-04-26/testrun_{0}.npz'.format(outname))
            A = f['A_best']
            s = np.max(f['s_history'])

            d_preds, est = predictive_bayesian_bdeu(A, 5, pred_data, tr_data)

            norm_c = np.sum(d_preds)
            #print(norm_c)

            d_preds_n  = d_preds/norm_c
            #print(np.sum(d_preds_n))

            scale = 1.0
            if USE_SCALES:
                scale = np.abs(sbest) / np.abs(s)
            sum_scales += scale
            for pi, p in enumerate(d_preds_n):
                #print(p)
                all_preds[pi] += p * scale
                pass
            #print(np.max(d_preds_n))


            f.close()

    # all_preds = all_preds / sum_scales # unnecessary because we normalize anyway

    all_norm = np.sum(all_preds)
    for ap in all_preds:
        print(ap/all_norm)


def process_disco():

    keys, tr_data = load_data()

    with open('{0}/data/test.txt'.format(curdir)) as f:
        #
        data_raw = [ l.split() for l in f]
    pred_data = [ [int(d) for d in l] for l in data_raw[1:] ]
    pred_data = np.array(pred_data)

    f = np.load('../disco-output/out-long-mi7.txt_04.npz')
    A = f['A']
    f.close()


    d_preds, est = predictive_bayesian_bdeu(A, 20, pred_data, tr_data)

    all_preds = np.zeros(pred_data.shape[0])
    for pi, p in enumerate(d_preds):
        #print(p)
        all_preds[pi] = p

    all_norm = np.sum(all_preds)
    for ap in all_preds:
        print(ap/all_norm)


if __name__ == "__main__":
    #
    process_python_ensembles()
    #process_disco()

    #f = np.load('../disco-output/out-long-mi7.txt_04.npz')
    #A = f['A']
    #f.close()
    #print(A)

    #stuff.visualize_network(A, fname='mcmc04')


