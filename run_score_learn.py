#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Project in probabilistic models
# Aaro Salosensaari

# runner for score learner

import numpy as np
import scipy.special

import warnings

np.random.seed(42)
np.seterr(all='warn')


curdir = '..'

from score_structure_learn import run_bic_climber


if __name__ == "__main__":
    for outletter in "ABCDEFGH":
        outname = "bic_dnp_2500_{0}".format(outletter)
        A_best, s_best, s_history, A_init, bic_valid_init, bic_valid, ll_valid_init, ll_valid = run_bic_climber(outname)
        print('best train score', s_best)
        print('valid bic scores: init, best', bic_valid_init, bic_valid)
        print('valid ll scores: init, best', ll_valid_init, ll_valid)
        np.savez(curdir+'/testoutput-2017-04-12b/testrun{0}'.format(outname), A_best=A_best, s_history=s_history, A_init=A_init, bic_valid=bic_valid, ll_valid=ll_valid)
    pass

