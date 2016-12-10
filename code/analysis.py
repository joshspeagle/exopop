#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script contains functions that run the analysis
on an input catalog (real or simulated).

All of the MCMC tuning parameters are hard coded. Sorry!
"""

from __future__ import division, print_function

import os
import h5py
import corner
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as pl

from load_data import load_candidates, load_detection_efficiency
from population import ProbabilisticModel, Dataset, Population


def inverse_detection_efficiency(pop, censor, catalog, truth=None):

    c = np.log(catalog) # log-ify catalog
    ind = [np.digitize(x, b) for x, b in zip(c.T, censor.bins)] # grab bin locations for objects on efficiency map
    weights = np.exp(-censor.lnprob[ind]) # compute weights
    val, x, y = np.histogram2d(c[:, 0], c[:, 1], pop.bins, weights=weights) # bin with raw weight (mean map)
    var, x, y = np.histogram2d(c[:, 0], c[:, 1], pop.bins, weights=weights**2) # bin with weight**2 (variance map)

    val[~np.isfinite(val)] = 0.0 # failures = 0
    var[~np.isfinite(var)] = 0.0 # failures = 0

    # build the model
    v = pop.initial() # initialize values
    lg = np.log(val).flatten() # log(occurence rate)
    m = np.isfinite(lg) # mask
    lg[~m] = 0.0 # set to 0 for failures
    v[-len(lg):] = lg # reorder population

    # compute marginalized errors
    marg = [np.sum(val, axis=(i+1) % 2) for i in range(2)] # x,y marginalization
    norm = [np.sum(a * np.diff(pop.bins[i])) for i, a in enumerate(marg)] # normalization
    literature = [(pop.bins[i], a/norm[i],
                   np.sqrt(np.sum(var, axis=(i+1) % 2))/norm[i])
                  for i, a in enumerate(marg)] # normalized histograms

    # extrapolate to Earth
    o = np.argsort(c[:, 0]) # order periods
    s = c[o, :] # sorted data (by period)
    ws = weights[o] # sorted weights
    m = np.isfinite(ws) * (s[:, 1] <= np.log(2)) * (s[:, 1] >= np.log(1)) # mask for selection box (1<M/M_Earth<2)
    cs = np.cumsum(ws[m]) # cdf
    cs_var = np.cumsum(ws[m] ** 2) # cdf(error)
    i = s[m, 0] > np.log(50) # limit to periods above 50 days

    # fit a line
    A = np.vander(s[m, 0][i], 2) #  Vandermonde matrix (x^1, x^2)
    Cinv = np.diag(1.0 / cs_var[i]) # precision matrix (C^-1)
    S = np.linalg.inv(np.dot(A.T, np.dot(Cinv, A))) # conditioned precision matrix (multivariate normal)
    mu = np.dot(S, np.dot(A.T, np.dot(Cinv, cs[i]))) # mean vector

    return v, val, var, literature, (mu, S), (s[m, 0], cs, cs_var)
