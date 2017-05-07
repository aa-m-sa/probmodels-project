#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Project in probabilistic models
# Aaro Salosensaari

# runner for score learner

import numpy as np
import scipy.special

import warnings
from datetime import datetime

#np.random.seed(4267)
#np.random.seed(4269)
np.random.seed(9999)
np.seterr(all='warn')


curdir = '..'

from score_structure_learn import run_bdeu_climber


if __name__ == "__main__":
    for outletter in "IJKLMNOPQRST":
        print(outletter)
        print(str(datetime.now()))
        outname = "bdeu_newfix_2500_{0}".format(outletter)
        A_best, s_best, s_history, A_init  = run_bdeu_climber(outname)
        print('best train score', s_best)
        np.savez(curdir+'/testoutput-2017-04-26/testrun_{0}'.format(outname), A_best=A_best, s_history=s_history, A_init=A_init)
    pass


