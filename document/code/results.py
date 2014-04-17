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
from matplotlib.ticker import FormatStrFormatter

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

samples = samples[-200000:, :]  # [::50, :]

# Load the true gamma earth if it exists.
fn = os.path.join(bp, "gamma.txt")
if os.path.exists(fn):
    gamma = [np.exp(float(open(fn).read())) / 42557.0]
else:
    gamma = None  # [5.7 / 100, 1.7 / 100, 2.2 / 100]

# Load the extrapolated value.
ext = np.array(open(os.path.join(bp, "extrap.txt"), "r").read().split(),
               dtype=float) / 42557.0
ext /= (np.log(400) - np.log(200)) * (np.log(2) - np.log(1))

# Compute and plot gamma_earth.
rates = pop.get_lnrate(samples, [np.log(365.), np.log(1.0)])
fracs = rates - np.log(42557.0)
a, b, c = triangle.quantile(fracs, [0.16, 0.5, 0.84])
print("{0}^{{+{1}}}_{{-{2}}}".format(b, c-b, b-a))

fig = pl.figure()
ax = fig.add_subplot(111)
ax.hist(fracs, 50, color="k", histtype="step", normed=True)

if gamma is not None:
    ax.axvline(np.log(gamma[0]), color="k", alpha=0.5, lw=3)

ax.axvline(np.log(ext[0]+ext[1]), color="k", ls="dashed")
ax.axvline(np.log(ext[0]+ext[2]), color="k", ls="dashed")
ax.axvline(np.log(ext[0]), color="k")

# ax.axvline(b, color="k")
# ax.axvline(c, color="k", ls="dashed")
# ax.axvline(a, color="k", ls="dashed")

ax.set_xlabel(r"$\ln \Gamma_\oplus$")
ax.set_ylabel(r"$p(\ln \Gamma_\oplus)$")

a2 = ax.twiny()
a2.set_xlim(100 * np.exp(ax.get_xlim()))
a2.set_xscale("log")
# a2.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
a2.set_xlabel(r"$\Gamma_\oplus\,[\%]$")

fig.savefig(os.path.join(bp, "rate.png"))
fig.savefig(os.path.join(bp, "rate.pdf"))

# Plot some posterior samples of the rate function.
somesamples = samples[np.random.randint(len(samples), size=50), :]
fig = pop.plot_2d(somesamples, censor=model.censor, catalog=np.log(catalog),
                  err=err, true=truth, labels=labels, top_axes=top_axes,
                  literature=literature)
fig.savefig(os.path.join(bp, "results.png"))
fig.savefig(os.path.join(bp, "results.pdf"))
