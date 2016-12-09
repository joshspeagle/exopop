#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script runs the analysis on either the real catalog or a simulated one.
To run on the real data, call:

    python main.py

and to run on a simulation, call:

    python main.py path/to/the/simulated/catalog/

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


def inverse_detection_efficiency(pop, censor, catalog, err, truth=None):
    c = np.log(catalog)
    ind = [np.digitize(x, b) for x, b in zip(c.T, censor.bins)]
    weights = np.exp(-censor.lnprob[ind])
    val, x, y = np.histogram2d(c[:, 0], c[:, 1], pop.bins, weights=weights)
    var, x, y = np.histogram2d(c[:, 0], c[:, 1], pop.bins, weights=weights**2)

    val[~np.isfinite(val)] = 0.0
    var[~np.isfinite(var)] = 0.0

    # Build the model for plotting.
    v = pop.initial()
    lg = np.log(val).flatten()
    m = np.isfinite(lg)
    lg[~m] = 0.0
    v[-len(lg):] = lg

    # Compute the marginalized errorbars.
    marg = [np.sum(val, axis=(i+1) % 2) for i in range(2)]
    norm = [np.sum(a * np.diff(pop.bins[i])) for i, a in enumerate(marg)]
    literature = [(pop.bins[i], a/norm[i],
                   np.sqrt(np.sum(var, axis=(i+1) % 2))/norm[i])
                  for i, a in enumerate(marg)]

    # Plot the results.
    labels = ["$\ln T/\mathrm{day}$", "$\ln R/R_\oplus$"]
    top_axes = ["$T\,[\mathrm{days}]$", "$R\,[R_\oplus]$"]
    fig = pop.plot_2d(v, censor=censor, catalog=np.log(catalog),
                      err=err, true=truth, labels=labels,
                      top_axes=top_axes, literature=literature)

    # Extrapolate to the Earth.
    o = np.argsort(c[:, 0])
    s = c[o, :]
    ws = weights[o]
    m = np.isfinite(ws) * (s[:, 1] <= np.log(2)) * (s[:, 1] >= np.log(1))
    cs = np.cumsum(ws[m])
    cs_var = np.cumsum(ws[m] ** 2)
    i = s[m, 0] > np.log(50)

    # Do the linear fit.
    A = np.vander(s[m, 0][i], 2)
    Cinv = np.diag(1.0 / cs_var[i])
    S = np.linalg.inv(np.dot(A.T, np.dot(Cinv, A)))
    mu = np.dot(S, np.dot(A.T, np.dot(Cinv, cs[i])))

    # Compute the predictive value.
    ys = np.dot(np.array([[np.log(200), 1],
                          [np.log(400), 1]]),
                np.random.multivariate_normal(mu, S, 5000).T)
    frac = np.diff(ys, axis=0)
    q = corner.quantile(frac, [0.16, 0.5, 0.84])

    fig2 = pl.figure()
    ax = fig2.add_subplot(111)
    a = np.vander(np.linspace(np.log(50), np.log(400), 500), 2)
    y = np.dot(a, np.random.multivariate_normal(mu, S, 50).T)
    ax.plot(a[:, 0], y, "r", alpha=0.3)
    ax.plot(a[:, 0], np.dot(a, mu), "--r", lw=2)
    ax.errorbar(s[m, 0], cs, yerr=np.sqrt(cs_var), fmt="k", capsize=0)

    return val, var, literature, fig, fig2, (q[1], np.diff(q))


def main(bp, real_data, ep_bins=False, ignore_uncert=False):
    if ignore_uncert:
        bp = "{0}-no-uncert".format(bp)
    try:
        os.makedirs(bp)
    except os.error:
        pass

    # Load the censoring function.
    censor = load_detection_efficiency()

    # The values from EP's paper (+some made up numbers).
    lpb, lrb = censor.bins
    x, y = lpb[::4], lrb[::4]

    # Load the candidates.
    if real_data:
        ids, catalog, err = load_candidates()
        truth = None
        m = np.log(catalog[:, 0]) > np.min(x)
        catalog, err = catalog[m], err[m]
        if ignore_uncert:
            err = np.zeros_like(err)
    else:
        catalog, err, truth = \
            pickle.load(open(os.path.join(bp, "catalog.pkl")))
    K = 1 if ignore_uncert else 512
    dataset = Dataset.sample(catalog, err, samples=K, censor=censor,
                             functions=[np.log, np.log])
    print("{0} entries in catalog".format(dataset.catalogs.shape[1]))

    # ...
    rerr = [np.log(catalog[:, 1]) - np.log(catalog[:, 1]-err[:, 1]),
            np.log(catalog[:, 1]+err[:, 1]) - np.log(catalog[:, 1])]
    err = [0, rerr]

    # Build the binned model.
    bins = [x, y]
    print("Run inference on a grid with shape: {0}"
          .format([len(b)-1 for b in bins]))
    pop = Population(bins, censor.bins, 11.0)
    model = ProbabilisticModel(dataset, pop, censor, [3.6, 2.6, 1.6, 0.0],
                               np.array([2.0, 0.5, 0.3, 0.3]) / 2.4)

    # Do V-max.
    val, var, literature, fig1, fig2, ext = \
        inverse_detection_efficiency(pop, censor, catalog, err, truth)

    open(os.path.join(bp, "extrap.txt"), "w").write(
        "{0} -{1} +{2}".format(ext[0], *(ext[1])))

    # Plot the vmax results.
    labels = ["$\ln T/\mathrm{day}$", "$\ln R/R_\oplus$"]
    top_axes = ["$T\,[\mathrm{days}]$", "$R\,[R_\oplus]$"]
    fig1.savefig(os.path.join(bp, "vmax.png"))
    fig1.savefig(os.path.join(bp, "vmax.pdf"))

    # Plot the extrapolation.
    fig2.savefig(os.path.join(bp, "extrapolation.png"))
    fig2.savefig(os.path.join(bp, "extrapolation.pdf"))

    # Save the model and the other things needed for plotting the results.
    pickle.dump((model, catalog, err, truth, labels, top_axes,
                 literature),
                open(os.path.join(bp, "model.pkl"), "w"), -1)

    # Set up the output files.
    nblock = 500
    N, ndim, nhyper = 2000 * nblock, len(pop), 4

    samples = np.empty((nblock, ndim))
    hyper = np.empty((nblock, nhyper))
    lnprob = np.empty(nblock)

    fn = os.path.join(bp, "results.h5")
    with h5py.File(fn, "w") as f:
        f.create_dataset("samples", shape=(N, ndim), dtype=np.float64)
        f.create_dataset("hyper", shape=(N, nhyper), dtype=np.float64)
        f.create_dataset("lnprob", shape=(N,), dtype=np.float64)

    for i, (th, hy, lp, acc) in enumerate(model.sample()):
        n = i % nblock
        samples[n, :] = th
        hyper[n, :] = hy
        lnprob[n] = lp
        if n == nblock - 1:
            print(i+1, (i+1) / N, np.max(lnprob), acc)
            s = slice(i-n, i+1)
            with h5py.File(fn, "a") as f:
                f.attrs["iteration"] = i+1
                f["samples"][s, :] = samples
                f["hyper"][s, :] = hyper
                f["lnprob"][s] = lnprob

        if i >= N-1:
            break

    pl.clf()
    pl.plot(samples[:, 0])
    pl.savefig("test.png")

    pl.clf()
    pl.plot(lnprob)
    pl.savefig("test-lp.png")

    pl.clf()
    pl.plot(hyper[:, 0])
    pl.savefig("test-hyper.png")


if __name__ == "__main__":
    import sys
    ignore = False
    if "--ignore" in sys.argv:
        sys.argv.remove("--ignore")
        ignore = True
    if len(sys.argv) > 1:
        main(sys.argv[1], False, ignore_uncert=ignore)
    else:
        main("results", True, ignore_uncert=ignore)
