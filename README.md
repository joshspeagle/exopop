Exoplanet population inference
==============================

This repository contains the a forked version of the
code and text for the paper [Exoplanet population
inference and the abundance of Earth analogs from noisy, incomplete catalogs](
http://arxiv.org/abs/1406.3020)
by Daniel Foreman-Mackey, David W. Hogg, and Timothy D. Morton and submitted
to ApJ. Most of the modifications to the original code were done for a final
project for Harvard's Exoplanets astronomy course (ASTRO189) to make things
more amenable to running iteractively in an IPython Notebook.

The code lives in the `code` directory and the LaTeX source code for the
original paper is in `document`.

**Code**

The meat of the original probabilistic model is implemented in `population.py`. Then
there are a set of scripts available that can generate the figures from the
paper. Details are located in the docstrings but essentially:

* `simulate.py` generates synthetic catalogs from known occurrence rate
  density functions (`simulate_cat.py` is the interactive version),
* `main.py` does the MCMC analysis on either real or simulated catalogs (`analysis.py` is
  the interactive version), and
* `results.py` analyzes the results of the MCMC, makes some figures, and thins
  the chain to the published form (`[TBD].py` is the interactive version).
* The entire analysis framework is described and outlined in `EXOPOP.ipynb`, which contains
  all of the relevant figures used in my final project.

Results
-------

The original simulated catalogs and results are [available online on figshare](
http://dx.doi.org/10.6084/m9.figshare.1051864). The new catalogs used for my project
can be found in the `data\` folder.

Attribution
-----------

The original code is associated with and written specifically for [Foreman-Mackey, Hogg, & Morton (2014)](http://arxiv.org/abs/1406.3020).
If you make any use of it, please cite the original paper:

```
@article{exopop,
   author = {{Foreman-Mackey}, D. and {Hogg}, D.~W. and {Morton}, T.~D.},
    title = {Exoplanet population inference and the abundance of
             Earth analogs from noisy, incomplete catalogs},
  journal = {ArXiv --- submitted to ApJ},
     year = 2014,
   eprint = {1406.3020}
}
```

If you make any use of some of the new code contained here, please contact me.

License
-------

Copyright 2014 Daniel Foreman-Mackey (who pretty much did everything)
Copyright 2016 Joshua Speagle (who is immensely grateful for DFM providing this fantastic public code)

The code in this repository is made available under the terms of the MIT
License. For details, see the LICENSE file.
