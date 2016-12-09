#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This is a modification of simulate.py to interactively generate synthetic catalogs from known occurrence rate densities.

Run it with the "smooth" option to sample for the smooth model (Catalog A in
the paper) and without for Catalog B.

"""

from __future__ import division, print_function

try:
    from savefig import monkey_patch
except ImportError:
    pass
else:
    monkey_patch()

import os
import numpy as np
import cPickle as pickle
from itertools import product
from scipy.misc import logsumexp

from load_data import load_detection_efficiency, load_candidates
from population import SeparablePopulation, Histogram, BrokenPowerLaw




def sim(out_loc, fname, seed=None, state=None, smooth=False):
    
    if seed is not None:
        np.random.seed(seed) # intialize seed if not given
    try:
        os.makedirs(out_loc) # if directory doesn't exist, make it
    except os.error:
        pass

    if state is not None:
        np.random.set_state(state) # initialize state if not given

    # save random state (so we can reproduce things later)
    pickle.dump(np.random.get_state(),
                open(os.path.join(out_loc, fname+".state.pkl"), "wb"), -1)

    # load censoring function
    censor = load_detection_efficiency()

    # load candidates to get fractional uncertainties
    inds, mu, sig = load_candidates()
    ferr = sig[:, 1] / mu[:, 1]

    # values from Petigura+'s paper (plus some made up numbers from DFM)
    lpb, lrb = censor.bins
    x, y = lpb[::4], lrb[::4]
    p_vals = np.log(np.array([8.9, 13.7, 15.8, 15.2, 13.3, 12.2]))
    r_vals = np.log(np.array([11, 11.5, 12, 14.2, 18.6, 5.9, 1.9, 1, 0.9, 0.7,
                              0.5, 0.5]))

    # normalize underlying distributions
    p_vals -= logsumexp(p_vals + np.log(np.diff(x)))
    r_vals -= logsumexp(r_vals + np.log(np.diff(y)))

    # build a synthetic population
    norm = 10.5
    if smooth:
        truth = [norm, 0.5, -0.2, 4.0, 0.8, -1.5, 1.0]
        pdist = BrokenPowerLaw(lpb)
        rdist = BrokenPowerLaw(lrb)
    else:
        truth = np.concatenate([[norm], p_vals, r_vals])
        pdist = Histogram(x, base=lpb)
        rdist = Histogram(y, base=lrb)
    pop0 = SeparablePopulation([pdist, rdist], lnnorm=truth[0])
    open(os.path.join(out_loc, fname+".gamma.txt"), "w").write(
        "{0}".format(pop0.get_lnrate(truth, [np.log(365), 0.0])[0]))

    # sample from this censored population
    lnrate = np.array(censor.lnprob[1:-1, 1:-1])
    lnrate += pop0.evaluate(truth)
    catalog = np.empty((0, 2))
    for i, j in product(xrange(len(lpb)-1), xrange(len(lrb)-1)):
        area = (lpb[i+1] - lpb[i]) * (lrb[i+1] - lrb[i])
        k = np.random.poisson(np.exp(lnrate[i, j]) * area)
        if k == 0:
            continue
        entry = np.vstack((np.random.uniform(lpb[i], lpb[i+1], k),
                           np.random.uniform(lrb[j], lrb[j+1], k))).T
        catalog = np.concatenate((catalog, entry), axis=0)

    # add in observational uncertainties
    catalog = np.exp(catalog)
    i = np.random.randint(len(ferr), size=len(catalog))
    err = np.vstack([np.zeros(len(catalog)), ferr[i] * catalog[:, 1]]).T
    catalog += err * np.random.randn(*(err.shape))

    truth = [(pdist.base, np.exp(pdist(truth[1:1+len(pdist)]))),
             (rdist.base, np.exp(rdist(truth[1+len(pdist):])))]

    # save results
    pickle.dump((catalog, err, truth),
                open(os.path.join(out_loc, fname+".cat.pkl"), "wb"), -1)

    return len(catalog)
