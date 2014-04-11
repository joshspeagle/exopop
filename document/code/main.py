#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

try:
    from savefig import monkey_patch
except ImportError:
    pass
else:
    monkey_patch()

import os
import h5py
import emcee
import numpy as np
import cPickle as pickle
from itertools import product

from load_data import (transit_lnprob0, ln_period0, load_completenes_sim,
                       load_candidates)
from population import (CensoringFunction,
                        ProbabilisticModel, Dataset, Population,
                        SmoothPopulation)


def main(ep_bins=True):
    bp = "main"
    try:
        os.makedirs(bp)
    except os.error:
        pass

    # Bins chosen to include EP's bins.
    if ep_bins:
        per_rng, np_bins = np.log([6.25, 100.0]), 4*4
        rp_rng, nr_bins = np.log([1.0, 16]), 4*8
        p_vals = np.array([8.9, 13.7, 15.8, 15.2])
        r_vals = np.array([12, 14.2, 18.6, 5.9, 1.9, 1, 0.9, 0.7])
    else:
        per_rng, np_bins = np.log([6.25, 400.0]), 4*6
        rp_rng, nr_bins = np.log([0.5, 32]), 4*12
        p_vals = np.array([8.9, 13.7, 15.8, 15.2, 15., 14.8])
        r_vals = np.array([11, 11.5, 12, 14.2, 18.6, 5.9, 1.9, 1, 0.9, 0.7,
                           0.5, 0.5])

    # Load the data.
    ln_P_inj, ln_R_inj, recovered = load_completenes_sim(per_rng=per_rng,
                                                         rp_rng=rp_rng)
    tlp = lambda lnp, lnr: transit_lnprob0 - 2.*(lnp - ln_period0)/3

    # Set up the censoring function.
    censor = CensoringFunction(np.vstack((ln_P_inj, ln_R_inj)).T, recovered,
                               bins=(np_bins, nr_bins),
                               range=[per_rng, rp_rng],
                               transit_lnprob_function=tlp)

    # The values from EP's paper (+some made up numbers).
    lpb, lrb = censor.bins
    x, y = lpb[::4], lrb[::4]
    literature = [(x, p_vals / np.sum(p_vals*np.diff(x))),
                  (y, r_vals / np.sum(r_vals*np.diff(y)))]

    # Load the candidates.
    ids, catalog, err = load_candidates()
    dataset = Dataset.sample(catalog, err, samples=64, censor=censor,
                             functions=[np.log, np.log])
    print("{0} entries in catalog".format(len(catalog)))

    # Build the binned model.
    bins = [x, y]
    print("Run inference on a grid with shape: {0}"
          .format([len(b)-1 for b in bins]))
    pop = Population(bins, censor.bins, 11.0)
    pop = SmoothPopulation([7.0, 0.2, 0.0, 0.2], pop)
    model = ProbabilisticModel(dataset, pop, censor)

    # Compute the vmax histogram.
    ix0 = np.digitize(np.log(catalog[:, 0]), bins[0])
    iy0 = np.digitize(np.log(catalog[:, 1]), bins[1])
    ix = np.digitize(np.log(catalog[:, 0]), lpb)
    iy = np.digitize(np.log(catalog[:, 1]), lrb)
    lp = censor.lnprob[ix, iy]
    grid = np.zeros((len(bins[0])-1, len(bins[1])-1))
    counts = np.zeros_like(grid)
    var = np.zeros_like(grid)
    for i, j in product(range(len(bins[0])-1), range(len(bins[1])-1)):
        m = (ix0 == i+1) * (iy0 == j+1)
        counts[i, j] = np.sum(m)
        if counts[i, j] == 0:
            continue
        v = lp[m]
        grid[i, j] = np.sum(np.exp(-v[np.isfinite(v)]))
        var[i, j] = np.sum(np.exp(-2*v[np.isfinite(v)]))
    grid[np.isinf(grid)] = 0.0
    var[np.isinf(var)] = 0.0
    std = np.sqrt(var)

    # Turn the vmax numbers into something that we can plot.
    v = pop.initial()
    samples = np.zeros((24, len(v)))
    for i in range(samples.shape[0]):
        lg = np.log(grid+std*np.random.randn(*(var.shape))).flatten()
        m = np.isfinite(lg)
        lg[~m] = 0.0
        samples[i] = v
        samples[i, -len(lg):] = lg

    # Plot the vmax results.
    rerr = [np.log(catalog[:, 1]) - np.log(catalog[:, 1]-err[:, 1]),
            np.log(catalog[:, 1]+err[:, 1]) - np.log(catalog[:, 1])]
    labels = ["$\ln T/\mathrm{day}$", "$\ln R/R_\oplus$"]
    top_axes = ["$T\,[\mathrm{days}]$", "$R\,[R_\oplus]$"]
    fig = pop.plot_2d(samples, censor=censor, catalog=np.log(catalog),
                      err=[0, rerr],
                      labels=labels, top_axes=top_axes, literature=literature)
    fig.savefig(os.path.join(bp, "vmax.png"))

    # Save the model and the other things needed for plotting the results.
    pickle.dump((model, catalog, [0, rerr], labels, top_axes, literature),
                open(os.path.join(bp, "model.pkl"), "w"), -1)

    assert 0

    # Set up the sampler.
    p0 = pop.initial()
    print("Initial ln-prob = {0}".format(model.lnprob(p0)))
    ndim, nwalkers = len(p0), 200
    sampler = emcee.EnsembleSampler(nwalkers, ndim, model)

    # Initialize the walkers.
    pos = [p0 + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]
    print("Sampling {0} dimensions with {1} walkers".format(ndim, nwalkers))

    # Make sure that all the initial positions have finite probability.
    finite = np.isfinite(map(model.lnprob, pos))
    assert np.all(finite), "{0}".format(np.sum(finite))

    # Run the sampler.
    N = 50000
    fn = os.path.join(bp, "results.h5")
    with h5py.File(fn, "w") as f:
        f.create_dataset("chain", shape=(nwalkers, N, ndim), dtype=np.float64)
        f.create_dataset("lnprob", shape=(nwalkers, N), dtype=np.float64)
    for i, (p, lp, s) in enumerate(sampler.sample(pos, iterations=N)):
        if (i + 1) % 100 == 0:
            print(i+1, np.max(lp))
        with h5py.File(fn, "a") as f:
            f["chain"][:, i, :] = p
            f["lnprob"][:, i] = lp


if __name__ == "__main__":
    main(False)