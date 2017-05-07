#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Project in probabilistic models
# Aaro Salosensaari

# 2nd method for structure learning:
# order MCMC by Koller And Firedman 2003

import numpy as np
import scipy.special
import itertools

from datetime import datetime


from stuff import load_data

# network format:
# a np matrix A of N x N, A[i,j] = 1 -> arc from i to j
# i = corresp ith letter in alphabet A,B,C,...,Z

curdir = '..'

par_score_dyn_memory = dict() # log sum_{U in allowed parents} score (X_i, U | D)

par_local_dyn_memory = dict() # log; score(i, U | D)

itercache = dict()

eps = np.finfo(np.double).eps

print('setup...')
print(str(datetime.now()))
for f in range(27):
    fr = np.arange(f)
    for n in range(6):
        parents_inds = []
        for parents_tup in itertools.combinations(fr, n):
            parents_inds.append(np.array(parents_tup).astype(np.int))
        itercache[(f,n)] = parents_inds

print('setup done.')
print(str(datetime.now()))

def most_probable_parents_var(i, feasible_parents, max_indg):
    # we have calculated all necessary stuff during algorithm computation, just retrieve

    feasible_len = feasible_parents.shape[0]
    n = np.min((feasible_len, max_indg))

    s_max = -np.inf
    parents_max = None
    for n_parents in range(n+1):
        # sum over U_k = U_n_parents

        #ic = itertools.combinations(feasible_parents, n_parents)
        ic = itercache[(feasible_len, n_parents)]


        for parents_ind in ic:
            parents = feasible_parents[parents_ind]
            if (i, parents.tostring()) in par_local_dyn_memory:
                tmp  = par_local_dyn_memory[(i, parents.tostring())]
                if tmp > s_max:
                    s_max = tmp
                    parents_max = parents
            else:
                raise ValueError('{0}, {1} should be in par_local_dyn_memory!'.format(i, parents))
    return parents_max


def most_probable_parents_var_nocache(i, feasible_parents, data, max_indg, alpha):

    feasible_len = feasible_parents.shape[0]
    n = np.min((feasible_len, max_indg))
    data_len, n_vars = data.shape

    ri = 3 # constant for all i
    s_max = -np.inf
    parents_max = None
    for n_parents in range(n+1):
        # sum over U_k = U_n_parents

        #ic = itertools.combinations(feasible_parents, n_parents)
        ic = itercache[(feasible_len, n_parents)]


        for parents_ind in ic:
            parents = feasible_parents[parents_ind]
            qi_terms = np.power(3,n_parents)
            #tmp  = par_local_dyn_memory[(i, parents.tostring())]
            tmp = bdeu_local_score(i, data, data_len, parents, qi_terms, ri, alpha)
            if tmp > s_max:
                s_max = tmp
                parents_max = parents
    return parents_max

def most_probable_parents_order(order, max_indg):
    # fetch most probable parents
    parents = []
    for ivar, i in enumerate(order):
        p = most_probable_parents_var(ivar, order[:(i+1)], max_indg)
        parents.append(p)
    return p


def bdeu_score_over_feasible_pars(i, feasible_parents, data, max_indg, alpha):
    # exacty computation over all subsets of parents

    feasible_len = feasible_parents.shape[0]
    n = np.min((feasible_len, max_indg))

    ri = 3 # constant for all i

    data_len, n_vars = data.shape
    #alpha = np.double(alpha)

    global par_local_dyn_memory

    # itertools is terribad efficienty-wise, but easiest to use
    # TODO can be sped up dyn prog approach?
    #s = np.double(0.0)

    # a very silly way to calculate s_size
    s_size = 0

    for n_parents in range(n+1):
        # sum over U_k = U_n_parents

        qi_terms = np.power(3,n_parents)

        #ic = itertools.combinations(feasible_parents, n_parents)
        ic = itercache[(feasible_len, n_parents)]


        for parents_ind in ic:
            s_size += 1

    s_list = np.zeros(s_size, dtype=np.double)
    s_ind = 0
    s_max = -np.inf
    for n_parents in range(n+1):
        # sum over U_k = U_n_parents

        qi_terms = np.power(3,n_parents)

        #ic = itertools.combinations(feasible_parents, n_parents)
        ic = itercache[(feasible_len, n_parents)]


        for parents_ind in ic:
            parents = feasible_parents[parents_ind]
            if (i, parents.tostring()) in par_local_dyn_memory:
                # maybe we have evaluated this particular set earlier?
                #s += par_local_dyn_memory[(i, parents.tostring())]
                tmp  = par_local_dyn_memory[(i, parents.tostring())]
                if tmp > s_max:
                    s_max = tmp
                s_list[s_ind] = tmp
                s_ind += 1
                continue

            # if not:
            tmp = bdeu_local_score(i, data, data_len, parents, qi_terms, ri, alpha)


            if tmp > s_max:
                s_max = tmp
            s_list[s_ind] = tmp
            s_ind += 1
            par_local_dyn_memory[(i, parents.tostring())] = tmp
            #tmpe = np.exp(tmp)
            ##par_local_dyn_memory[(i, parents.tostring())] = tmpe
            #s += tmpe
    if not s_ind == s_size:
        print('warn sind')
    if not s_max > -np.inf:
        print('warn s_max neginf')

    s_arr = np.array(s_list) - s_max
    s_exp = np.exp(s_arr)
    s = np.sum(s_exp)
    s_log = np.log(s) + s_max
    return s_log


def bdeu_local_score(i, data, data_len, parents, qi_terms, ri, alpha):
    confs_seen = dict() # N_ij

    confs_seen_by_ival = dict() # N_ijk

    n_confs_seen = 0


    for n in range(data_len):
        ival = data[n,i]

        conf_orig = data[n,:][parents]
        conf = conf_orig.tostring() # 'hash'

        if conf not in confs_seen:
            confs_seen[conf] = 1
            n_confs_seen += 1
        else:
            confs_seen[conf] += 1

        if (conf, ival) not in confs_seen_by_ival:
            confs_seen_by_ival[(conf, ival)] = 1
        else:
            confs_seen_by_ival[(conf, ival)] += 1

    tmp = 0.0
    for conf in confs_seen:
        Nij = confs_seen[conf]
        alpha_ij = alpha/qi_terms
        alpha_ijk = alpha_ij/ri # ri constant for all i
        tmp_d = Nij + alpha_ij
        if tmp_d < 2*eps:
            print('warn, Nij+alpha_ij = 0')
            tmp_d = eps
        tmp_a = scipy.special.gammaln(alpha_ij) - scipy.special.gammaln(tmp_d)
        #tmp_a = scipy.special.gamma(alpha_ij) / scipy.special.gamma(tmp_d)
        tmp_b = 0.0
        for ival in range(3):
            if (conf, ival) not in confs_seen_by_ival:
                Nijk = eps
            else:
                Nijk = confs_seen_by_ival[(conf, ival)]
            Nijk = np.double(Nijk)

            tmp_gna_ga = scipy.special.gammaln(Nijk + alpha_ijk) - scipy.special.gammaln(alpha_ijk)
            #tmp_gna_ga = scipy.special.gamma(Nijk + alpha_ijk)/scipy.special.gamma(alpha_ijk)
            tmp_b += tmp_gna_ga

        tmp += tmp_a + tmp_b
    return tmp

def p_data_given_order(order, data, max_indg, alpha, first=False):
    """
    computes P(data | order)
    """
    # order: array of size N , listing order as [predecessor, successor] i.e. root first, leaves at the end
    # if node i is listed before node j: node i is a possible parent of j
    # max_indg: max in-degree (no of parents) per node
    # return log


    # TODO candidate set  / highest scoring families
    # TODO pruning

    data_len, N = data.shape

    magic = 0.0

    global par_score_dyn_memory
    s = np.double(0.0)
    i = 0
    for ivar in order:
        # s_tmp is in log
        if first:
            print(i, ivar)
            print(str(datetime.now()))
        feasible_parents = order[:(i+1)]
        if (ivar, feasible_parents.tostring()) in par_score_dyn_memory:
            s_tmp = par_score_dyn_memory[(ivar, feasible_parents.tostring())] # do we have evaluated this feasible set earlier?
        else:
            s_tmp =  bdeu_score_over_feasible_pars(ivar, feasible_parents, data, max_indg, alpha) + magic
            par_score_dyn_memory[(ivar, feasible_parents.tostring())] = s_tmp
            # TODO dyn prog trick Niinimäki thesis eq 2.6?
        s += s_tmp
        i += 1

    return s - magic*i


def mcmc_algorithm(data, alpha=5, max_indg=6, max_iter=10000, burnin=100, thin=100):
    """
    traditional metropolis-hastings over orders
    """
    print('start')

    data_len, n_vars = data.shape

    order = np.random.permutation(n_vars)

    n_samples = int(np.floor((max_iter - burnin) / thin) + 1)

    samples = np.zeros((n_samples, n_vars))
    p_history = np.zeros(max_iter+1)
    o_history = np.zeros((max_iter+1, n_vars))
    p_accept = np.zeros(max_iter+1)
    o_accept = np.zeros((max_iter+1, n_vars))

    i = 0
    si = 0
    p_old = p_data_given_order(order, data, max_indg, alpha, first=True)
    p_history[i] = p_old
    o_history[i,:] = order

    print('init calcd', p_old)
    print(order)

    while i < max_iter:
        # only swaps
        swp_i, swp_j = np.random.randint(0, n_vars, 2)
        order_new = order.copy()
        order_new[swp_i] = order[swp_j]
        order_new[swp_j] = order[swp_i]
        # TODO one can use the substraction trick to speed up computations
        # TODO also pruning, approximate methods?

        p_new = p_data_given_order(order_new, data, max_indg, alpha)

        #d = np.log(p_new) - np.log(p_old)
        d = p_new - p_old


        if d >= 0:
            order = order_new
            p_old = p_new
        elif np.random.rand() < np.exp(d):
            order = order_new
            p_old = p_new

        o_history[i+1,:] = order_new
        p_history[i+1] = p_new

        o_accept[i+1,:] = order
        p_accept[i+1] = p_old


        if i % 5 == 0:
            print(i)
            print(str(datetime.now()))
            print(p_new, p_old)
            print(order_new)
            print(order)

        if (i - burnin) % thin == 0:
            print("sampled")
            samples[si,:] = order
            si += 1



        i += 1

    samples[-1,:] = order

    return samples, p_history, o_history, p_accept, o_accept




def run_order_mcmc_learner(outname, ensemble_size=5, max_indg=7):
    keys, data = load_data()
    print('ensemble start')
    print(str(datetime.now()))
    for walker in range(ensemble_size):
        print('walker', walker, 'of', ensemble_size)
        d = datetime.now()
        print(str(d))
        samples, p_history, o_history, p_accept, o_accept = mcmc_algorithm(
            data, alpha=5, max_indg=max_indg, max_iter=1000,burnin=300, thin=50)
        print('walker done, constructing network...')
        print(str(datetime.now()))
        # network of most probable parents
        A = np.zeros((data.shape[1], data.shape[1]))
        forder = samples[-1,:]
        parents = most_probable_parents_order(forder, max_indg)
        for i, ivar in enumerate(forder):
            for p in parents[i]:
                A[p, ivar] = 1
        print('done, saving...')
        print(str(datetime.now()))
        np.savez(
            curdir+'/testoutput-order-mcmc/run_{0}_{1}_w{2}'.format(outname, str(d.date()), walker),
            samples=samples, p_history=p_history, o_history=o_history,
            p_accept=p_accept, o_accept=o_accept, A_final = A)



    np.savez(curdir+'/testoutput-order-mcmc/dyn_mem_{0}_{1}'.format(outname, str(d.date())),
             par_score_dyn_memory=par_score_dyn_memory,
             par_local_dyn_memory=par_local_dyn_memory)
    # stuff

def test():
    keys, data = load_data()
    print(str(datetime.now()))
    p1 = p_data_given_order(np.array([17, 20, 16,  3,  5, 21, 22, 19,  8, 10,  6, 23, 24, 11,  4, 18,  14,  9,  1,  7, 25, 15,  0,  2, 12, 13]), data, 3, 5)
    print(p1)
    print(str(datetime.now()))
    p2 = p_data_given_order(np.array([17, 20, 16,  3,  5, 21, 22, 19,  8, 10,  6, 23, 24, 11,  4, 18,  14,  9,  1,  7, 25, 15,  0,  2, 12, 13]), data, 3, 5)
    print(p2)
    print(str(datetime.now()))
    p3 = p_data_given_order(np.array([17, 20, 16,  11,  5, 21, 22, 19,  8, 10,  6, 23, 24, 3,  4, 18,  14,  9,  1,  7, 25, 15,  0,  2, 12, 13]), data, 3, 5)
    print(p3)
    print(str(datetime.now()))

def test_samples():
    keys, data = load_data()
    # from disco
    samples_rwa = [
        '5 10 2 25 17 3 22 18 1 16 20 12 19 15 0 21 13 4 6 14 24 7 9 8 23 11',
        '5 10 2 25 17 3 22 18 13 20 19 12 16 1 15 21 0 4 7 14 24 6 9 8 23 11',
        '25 5 2 10 17 3 22 18 1 20 19 15 12 16 21 13 0 7 14 24 4 6 9 8 23 11',
        '25 5 2 10 17 3 22 13 1 12 18 15 20 16 21 14 0 7 19 24 4 9 11 8 23 6',
        '5 25 2 17 10 3 22 18 1 15 19 12 20 16 21 13 0 7 14 24 4 9 11 8 23 6',
        '5 25 2 17 10 3 22 18 16 13 20 12 19 1 0 15 21 4 8 24 23 9 11 7 6 14',
        '5 25 2 17 1 3 22 18 16 15 10 12 19 13 20 0 21 4 8 24 23 9 11 7 6 14',
        '13 25 2 17 1 3 22 18 5 15 10 12 19 16 20 0 21 4 6 24 23 9 8 14 11 7',
        '13 25 2 17 12 3 22 18 15 5 10 1 19 16 20 0 21 4 6 23 24 9 7 14 11 8',
        '13 25 12 17 2 3 22 18 15 5 19 1 10 20 16 0 21 4 8 11 24 9 6 7 23 14',
        '13 25 1 17 2 3 22 18 16 12 19 5 10 21 15 0 20 4 14 11 24 9 23 7 6 8',
        '17 25 1 13 2 3 22 18 19 12 16 5 10 21 15 0 4 7 14 24 23 9 8 11 20 6' ]

    samples = [[int(i) for i in s.split()] for s in samples_rwa]
    s = samples[0]
    for ii, i in enumerate(s):
        print(ii)
        feasible_parents = np.array(s[:ii], dtype=int)
        p = most_probable_parents_var_nocache(i, feasible_parents, data, 5, 5)
        # still too slow, I'm probably missing something
    print(p)


if __name__ == "__main__":
    #test()
    #run_order_mcmc_learner('test1', max_indg=5)
    test_samples()
