#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from score_structure_learn import load_data
from datetime import datetime

import warnings

import stuff

curdir = '..'


if __name__ == "__main__":

    f = np.load('../testoutput-2017-04-26/testrun_bdeu_newfix_2500_T.npz')
    A = f['A_best']
    f.close()
    #print(A)

    stuff.visualize_network(A, fname='bdeu2500')
