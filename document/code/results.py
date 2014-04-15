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
import sys
import h5py
import triangle
import numpy as np
import cPickle as pickle
import matplotlib.pyplot as pl

import load_data

bp = sys.argv[1]
model, catalog, err, truth, labels, top_axes, literature = \
    pickle.load(open(os.path.join(bp, "model.pkl")))
pop = model.population

with h5py.File(os.path.join(bp, "results.h5")) as f:
    i = int(f.attrs["iteration"])
    samples = f["samples"][:i, :]
    hyper = f["hyper"][:i, :]
    lnprob = f["lnprob"][:i]

for i in range(hyper.shape[1]):
    pl.clf()
    pl.plot(hyper[:, i])
    pl.savefig(os.path.join(bp, "time-hyper-{0:03d}.png".format(i)))

samples = samples[-200000:, :][::50, :]

# Load the true gamma earth if it exists.
fn = os.path.join(bp, "gamma.txt")
if os.path.exists(fn):
    gamma = [np.exp(float(open(fn).read())) / 42557.0]
else:
    gamma = [5.7 / 100, 1.7 / 100, 2.2 / 100]

# Load the extrapolated value.
ext = np.array(open(os.path.join(bp, "extrap.txt"), "r").read().split(),
               dtype=float) / 42557.0

# Compute and plot gamma_earth.
rates = pop.get_lnrate(samples, [np.log(365.), np.log(1.0)])
fracs = rates - np.log(42557.0)
a, b, c = triangle.quantile(fracs, [0.16, 0.5, 0.84])
print("{0}^{{+{1}}}_{{-{2}}}".format(b, c-b, b-a))
pl.clf()
pl.hist(fracs, 50, color="k", histtype="step", normed=True)

pl.gca().axvline(np.log(gamma[0]), color="r")
if len(gamma) > 1:
    pl.gca().axvline(np.log(gamma[0]+gamma[1]), color="r", ls="dashed")
    pl.gca().axvline(np.log(gamma[0]-gamma[2]), color="r", ls="dashed")

pl.gca().axvline(np.log(ext[0]+ext[1]), color="b", ls="dashed")
pl.gca().axvline(np.log(ext[0]+ext[2]), color="b", ls="dashed")
pl.gca().axvline(np.log(ext[0]), color="b")

pl.gca().axvline(b, color="k")
pl.gca().axvline(c, color="k", ls="dashed")
pl.gca().axvline(a, color="k", ls="dashed")

pl.savefig(os.path.join(bp, "rate.png"))

# Plot some posterior samples of the rate function.
somesamples = samples[np.random.randint(len(samples), size=50), :]
fig = pop.plot_2d(somesamples, censor=model.censor, catalog=np.log(catalog),
                  err=err, true=truth, labels=labels, top_axes=top_axes,
                  literature=literature)
fig.savefig(os.path.join(bp, "results.png"))
# fig.savefig(os.path.join(bp, "results.pdf"))
