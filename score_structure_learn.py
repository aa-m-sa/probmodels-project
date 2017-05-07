#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Project in probabilistic models
# Aaro Salosensaari

# 1st attempt at structure learning

# optimize BIC/AIC/BDeu score with greedy hill climbing

import numpy as np
import scipy.special

import warnings

from stuff import load_data
from stuff import initial_network
from stuff import initial_network_indg
from stuff import initial_network_in_outdg
from stuff import check_dagness

#np.random.seed(42)
np.seterr(all='warn')


curdir = '..'

# network format:
# a np matrix A of N x N, A[i,j] = 1 -> arc from i to j
# i = corresp ith letter in alphabet A,B,C,...,Z

ll_dyn_memory = dict()
bdeu_dyn_memory = dict()

def bdeu_score(A, data, alpha, dynamic=False, batch=False, dynamic_batch=None):
    """
    assumes no prior ie. log P(D|G)

    """
    n_parents = np.sum(A,0)

    qi_terms = np.power(3,n_parents)
    ri = 3 # constant for all i

    N = A.shape[0] # note, number of variables
    data_len = data.shape[0]
    alpha = np.double(alpha)

    s = np.double(0.0)
    for i in range(N):
        parents = np.where(A[:,i])[0]
        if  dynamic and (i, parents.tostring()) in bdeu_dyn_memory:
            if batch:
                tmp = bdeu_dyn_memory[(dynamic_batch, i, parents.tostring())]
            else:
                tmp = bdeu_dyn_memory[(i, parents.tostring())]
        else:
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

            tmp = np.double(0.0)
            for conf in confs_seen:
                Nij = confs_seen[conf]
                alpha_ij = alpha/qi_terms[i]
                alpha_ijk = alpha_ij/ri # ri constant for all i
                tmp_d = Nij + alpha_ij
                if tmp_d < np.finfo(np.double).eps:
                    print('warn, Nij+alpha_ij = 0')
                    tmp_d = np.finfo(np.double).eps
                tmp_a = scipy.special.gammaln(alpha_ij) - scipy.special.gammaln(tmp_d)
                tmp_b = np.double(0.0)
                for ival in range(3):
                    if (conf, ival) not in confs_seen_by_ival:
                        Nijk = np.finfo(np.double).eps
                    else:
                        Nijk = confs_seen_by_ival[(conf, ival)]
                    Nijk = np.double(Nijk)

                    #tmp_gna_ga = np.log(scipy.special.gamma(Nijk + alpha_ijk)/scipy.special.gamma(alpha_ijk))
                    tmp_gna_ga = scipy.special.gammaln(Nijk + alpha_ijk) - scipy.special.gammaln(alpha_ijk)
                    tmp_b += tmp_gna_ga
                tmp += tmp_a + tmp_b
                #with warnings.catch_warnings():
                #    warnings.filterwarnings('error')
                #    try:
                #        tmp += tmp_a + tmp_b
                #    except Warning:
                #        print('err')
                #        print(tmp_a)
                #        print(tmp_b)
                #        print(alpha,alpha_ij,alpha_ijk,Nij+alpha_ij, tmp_c, tmp_a)
                #        raise
            if dynamic:
                if batch:
                    bdeu_dyn_memory[(dynamic_batch, i, parents.tostring())] = tmp
                else:
                    bdeu_dyn_memory[(i, parents.tostring())] = tmp
        s += tmp
    return s #TODO return s => maximizers found empty network always. maybe we should negate it?
    # however no idea if it's theoretically justified in any way, probably something wrong in comp
    # maybe we are just sensitive to alpha (too small alpha, network ege count small, see lect slide p.24)
    # it appears that alpha = 10.0**20 might even work?

def ll_score(A, data, dynamic=False, batch=False, dynamic_batch=None):
    # LL term
    # sum_{i=1}^n sum_{j=1}^{q_i} sum_{k=1}^{r_i} N_{ijk} log(N_{ijk}/N_{ij}
    # notice that this is a sum of local scores
    # -> each change during a search is only a couple of steps at a time
    # -> most of the scores could be stored and retrieved on the next calc round?
    # note that ^will not 'sync' with batch based learning (data set will be different on each call)

    # number of parents for each node i

    n_parents = np.sum(A,0)

    N = A.shape[0]
    data_len = data.shape[0]

    ll = 0.0
    # replicate this with C/Cython? arrays?
    for i in range(N):
        parents = np.where(A[:,i])[0]
        if  dynamic and (i, parents.tostring()) in ll_dyn_memory:
            if batch:
                tmp = ll_dyn_memory[(dynamic_batch, i, parents.tostring())]
            else:
                tmp = ll_dyn_memory[(i, parents.tostring())]
        else:
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
            with warnings.catch_warnings():
                warnings.filterwarnings('error')
                for conf in confs_seen:
                    Nij = confs_seen[conf]
                    for ival in range(3):
                        if (conf, ival) not in confs_seen_by_ival:
                            Nijk = np.finfo(np.double).eps
                            #print('Warning: BIC score had to resort to eps for Nijk')
                            #print(parents)
                            #print(np.fromstring(conf, dtype=int))
                            #print(ll, i, ival)
                        else:
                            Nijk = confs_seen_by_ival[(conf, ival)]
                        Nijk = np.double(Nijk)
                        try:
                            tmp += Nijk * np.log(Nijk/Nij)
                        except Warning:
                            print('RuntimeWarning in BIC score')
                            print(Nijk, Nij)
                            print(ll, i)
            if dynamic:
                if batch:
                    ll_dyn_memory[(dynamic_batch, i, parents.tostring())] = tmp
                else:
                    ll_dyn_memory[(i, parents.tostring())] = tmp
        ll += tmp
    return ll

def bic_score(A, data, dynamic=False, batch=False, dynamic_batch=False):
    """
    BIC score for network structure A given data.

    Numerical evaluation not optimized much.
    # main loops are very parallizable
    # consider some kind of parallel solution + C/Fortran?
    """
    N = A.shape[0]
    data_len = data.shape[0]

    # denote
    # r_i = number of states of variable i == 3 always
    # q_i = number of possible configurations of parents = 3^{number of parents of i in graph A}

    # N_{ijk} = number of instances in data where variable i has value k and parents are in configuration j
    # N_{ij} = previous summed over k = 1 .. r_i


    ll = ll_score(A, data, dynamic=dynamic, batch=batch, dynamic_batch=dynamic_batch)


    # penalty term

    # sum_{i=1}^n (r_i -1)q_i
    n_parents = np.sum(A,0)

    # number of states its parents can be for each i
    qi_terms = np.power(3,n_parents)

    complexity = 2*np.sum(qi_terms)
    penalty = 0.5 * np.log(N) * complexity

    return ll - penalty


def aic_score(A, data, dynamic=False, batch=False, dynamic_batch=False):
    N = A.shape[0]
    data_len = data.shape[0]

    # denote
    # r_i = number of states of variable i == 3 always
    # q_i = number of possible configurations of parents = 3^{number of parents of i in graph A}

    # N_{ijk} = number of instances in data where variable i has value k and parents are in configuration j
    # N_{ij} = previous summed over k = 1 .. r_i


    ll = ll_score(A, data, dynamic=dynamic, batch=batch, dynamic_batch=dynamic_batch)


    # penalty term

    # sum_{i=1}^n (r_i -1)q_i
    n_parents = np.sum(A,0)

    # number of states its parents can be for each i
    qi_terms = np.power(3,n_parents)

    complexity = 2*np.sum(qi_terms)
    penalty = complexity

    return ll - penalty


def random_perturb(A_orig, perturb=20):
    """
    perturb network A, with size perturb

    don't care about in/outdg constraints, perturb assumed small anyway
    """
    n_add  = np.random.randint(perturb)
    n_remove = perturb - n_add
    # this ns are max to remove instead of precise amount removed

    # first remove n_remove random edges
    A = A_orig.copy()
    N = A.shape[0]

    n_r = 0
    for i in np.random.permutation(N):
        for j in np.random.permutation(N):
            if not A[i,j]:
                continue
            if n_r > n_remove:
                continue
            A[i,j] = 0
            n_r += 1

    # add n_add random edges, ensuring DAGness
    n_a = 0
    attempts = 0
    while n_a < n_add or attempts < 50:
        for i in np.random.permutation(N):
            for j in np.random.permutation(N):
                if A[i,j]:
                    continue
                if n_a > n_add:
                    continue
                A_try = A.copy()
                A_try[i,j] = 1
                if not check_dagness(A_try, i,j):
                    continue
                A = A_try
                n_a += 1
        attempts += 1

    return A




def find_nb_max(A_start, s_start, score, max_indg, max_outdg, tabu):
    # try deleting edges
    N = A_start.shape[0]
    A = A_start
    s = s_start
    # find best network starting from A by deleting edges
    A_dbest = None
    s_dbest = -np.inf
    for i in range(N):
        for j in range(N):
            if not A[i,j]:
                continue
            A_try = A.copy()
            A_try[i,j] = 0
            s_try = score(A_try)
            if  s_try > s_dbest:
                if A_try.tostring() in tabu:
                    continue
                s_dbest = s_try
                A_dbest = A_try

    # find best network starting from A by adding edges
    A_abest = None
    s_abest = -np.inf
    for j in range(N):
        if np.sum(A[:,j]) >= max_indg:
            continue
        for i in range(N):
            if A[i,j]:
                continue
            if np.sum(A[i,:]) >= max_outdg:
                continue
            A_try = A.copy()
            A_try[i,j] = 1
            # check that networks stays DAG
            if not check_dagness(A_try, i,j):
                continue
            s_try = score(A_try)
            if  s_try > s_abest:
                if A_try.tostring() in tabu:
                    continue
                s_abest = s_try
                A_abest = A_try

    # TODO swapping directions

    if s_abest > s_dbest:
        return A_abest, s_abest
    elif s_dbest > s_start:
        return A_dbest, s_dbest
    else:
        return A_start, s_start


def greedy_hill_climber(score, max_indg=5, max_outdg=5, max_iter=1000, outname=""):
    """
    score: larger is better; function accepts only one param (network)
    uses totally random restarts
    """
    # generate initial random network
    A_init = initial_network_indg(max_indg)

    A_prev = A_init
    s_prev = score(A_prev)

    A_best = A_prev
    s_best = s_prev

    tabu = set()

    s_history = np.zeros(max_iter)

    i = 0
    while i < max_iter:
        A_local_max, s_local_max = find_nb_max(A_prev, s_prev, score, max_indg, max_outdg, tabu)
        if s_prev >= s_local_max:
            A = initial_network_indg(max_indg)
            s = score(A)
        else:
            A = A_local_max
            s = s_local_max

        if s > s_best:
            A_best = A
            s_best = s

        tabu.add(A_prev.tostring())
        A_prev = A
        s_prev = s
        s_history[i] = s
        if i % 100 == 0:
            print(i, s_best, s)
            np.savez(curdir+"/testoutput/cur{0}_{1}".format(outname, str(i)), A)
            np.savez(curdir+"/testoutput/best{0}_{1}".format(outname, str(i)), A_best)

        i = i+1

    return A_best, s_best, s_history, A_init


def greedy_hill_climber_perturb(score, max_indg=5, max_outdg=5, max_iter=1000, outname="", perturb=50):
    """
    score: larger is better; function accepts only one param (network)
    uses small perturbs.
    """
    # generate initial random network
    A_init = initial_network_in_outdg(max_indg, max_outdg)

    A_prev = A_init
    s_prev = score(A_prev)

    A_best = A_prev
    s_best = s_prev

    tabu = set()

    s_history = np.zeros(max_iter)

    i = 0
    while i < max_iter:
        A_local_max, s_local_max = find_nb_max(A_prev, s_prev, score, max_indg, max_outdg, tabu)
        if s_prev >= s_local_max:
            A = random_perturb(A_prev, perturb)
            s = score(A)
        else:
            A = A_local_max
            s = s_local_max

        if s > s_best:
            A_best = A
            s_best = s

        tabu.add(A_prev.tostring())
        A_prev = A
        s_prev = s
        s_history[i] = s
        if i % 50 == 0:
            print(i, s_best, np.sum(A_best), s, np.sum(A))
            np.savez(curdir+"/testoutput/cur{0}_{1}".format(outname, str(i)), A)
            np.savez(curdir+"/testoutput/best{0}_{1}".format(outname, str(i)), A_best)

        i = i+1

    return A_best, s_best, s_history, A_init



def greedy_hill_climber_randombatch(score, data, max_indg=5, max_outdg=5, max_iter=1000, outname="", batch_size=100):
    """
    score: larger is better, accepts network,data,etc as params (for selecting data batch)
    uses totally random restarts

    consider only small batch of data at a time
    """
    # generate initial random network
    A_init = initial_network_indg(max_indg)

    N = data.shape[0]

    tmp = 0
    batches = []
    while tmp + batch_size < N:
        batches.append((tmp, tmp + batch_size))
        tmp += batch_size
    batches.append((tmp, N))

    n_batches = len(batches)


    A_prev = A_init
    s_prevs = np.zeros(n_batches)
    for bi,b in enumerate(batches):
        s_prevs[bi] = score(A_prev,data[b[0]:b[1]], dynamic=True, batch=True, dynamic_batch=bi)

    A_best = A_prev
    s_bests = s_prevs.copy()

    tabu = set()

    s_history = np.zeros((max_iter, n_batches))


    i = 0

    while i < max_iter:
        bi = np.random.randint(n_batches)
        b_start = batches[bi][0]
        b_end = batches[bi][1]

        score_batch = lambda x: score(x, data[b_start:b_end,:], dynamic=True, batch=True, dynamic_batch=bi)

        A_local_max, s_local_max = find_nb_max(A_prev, s_prevs[bi], score_batch, max_indg, max_outdg, tabu)

        if s_prevs[bi] >= s_local_max:
            A = initial_network_indg(max_indg)
            s = score_batch(A)
        else:
            A = A_local_max
            s = s_local_max

        if s > s_bests[bi]:
            A_best = A
            s_bests[bi] = s

        tabu.add(A_prev.tostring())
        A_prev = A
        s_prevs[bi] = s
        s_history[i,bi] = s
        if i % 10 == 0:
            print(i, s_bests, s)
            np.savez(curdir+"/testoutput/cur{0}_{1}".format(outname, str(i)), A)
            np.savez(curdir+"/testoutput/best{0}_{1}".format(outname, str(i)), A_best)

        i = i+1

    s_best = score(A_best, data)

    return A_best, s_best, s_history, A_init


def run_bic_climber(outname):
    """TODO: Docstring for main.
    :returns: TODO

    """
    global ll_dyn_memory
    ll_dyn_memory = dict()
    # todo consider saving / loading this?
    #
    keys, data = load_data()

    # first try: 2300 for tr, rest for valid
    #data_train = data[:2300]
    data_valid = data[2300:]
    data_train = data

    score = lambda x: bic_score(x, data_train, dynamic=True)
    #A_best, s_best, s_history, A_init = greedy_hill_climber(score, max_indg=3)
    #A_best, s_best, s_history, A_init = greedy_hill_climber(score, max_indg=5, max_outdg=5, max_iter=6000, outname=outname)
    A_best, s_best, s_history, A_init = greedy_hill_climber_perturb(score, max_indg=15, max_outdg=15, max_iter=500, outname=outname)

    #A_best, s_best, s_history, A_init = greedy_hill_climber_randombatch(bic_score, data, max_indg=5, outname=outname, batch_size=500)
    # batch is actually slower than full with ext memory...

    # evaluate network with data_valid
    bic_valid = bic_score(A_best, data_valid)
    bic_valid_init = bic_score(A_init, data_valid)
    ll_valid = ll_score(A_best, data_valid)
    ll_valid_init = ll_score(A_init, data_valid)

    return A_best, s_best, s_history, A_init, bic_valid_init, bic_valid, ll_valid_init, ll_valid

def run_aic_climber(outname):
    """TODO: Docstring for main.
    :returns: TODO

    """
    global ll_dyn_memory
    ll_dyn_memory = dict()
    # todo consider saving / loading this?
    #
    keys, data = load_data()

    # first try: 2300 for tr, rest for valid
    #data_train = data[:2300]
    data_valid = data[2300:]
    data_train = data

    score = lambda x: aic_score(x, data_train, dynamic=True)
    A_best, s_best, s_history, A_init = greedy_hill_climber_perturb(score, max_indg=15, max_outdg=15, max_iter=500, outname=outname)


    # evaluate network with data_valid
    aic_valid = aic_score(A_best, data_valid)
    aic_valid_init = aic_score(A_init, data_valid)
    ll_valid = ll_score(A_best, data_valid)
    ll_valid_init = ll_score(A_init, data_valid)

    return A_best, s_best, s_history, A_init, aic_valid_init, aic_valid, ll_valid_init, ll_valid


def run_bdeu_climber(outname):
    """
    note. still no CV
    """
    global ll_dyn_memory
    ll_dyn_memory = dict()
    # todo consider saving / loading this?
    #
    keys, data = load_data()

    # first try: 2300 for tr, rest for valid
    #data_train = data[:500]
    data_valid = data[2400:]
    data_train = data

    #alpha = 10.0**20
    #alpha = 10.**22
    alpha = 5.0
    print(np.log(alpha))
    #miter = 2000 # for full run 2000
    miter = 800
    score = lambda x: bdeu_score(x, data_train, alpha=alpha, dynamic=True)
    A_best, s_best, s_history, A_init = greedy_hill_climber_perturb(score, max_indg=26, max_outdg=26, max_iter=miter, outname=outname)


    # evaluate network with data_valid
    #aic_valid = aic_score(A_best, data_valid)
    #aic_valid_init = aic_score(A_init, data_valid)
    print(alpha)
    bd_valid = bdeu_score(A_best, data_valid, alpha=alpha)
    bd_valid_init = bdeu_score(A_init, data_valid, alpha=alpha)
    ll_valid = ll_score(A_best, data_valid)
    ll_valid_init = ll_score(A_init, data_valid,)
    print(bd_valid, bd_valid_init)
    print(ll_valid, ll_valid_init)

    return A_best, s_best, s_history, A_init
