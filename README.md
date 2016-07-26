MagniPy
==========

A repository for weak lensing magnification measurement code.

The notebook [measurement.ipynb](https://github.com/jesford/MagniPy/blob/master/measurement.ipynb) 
is a simple proof-of-concept analysis using the [TreeCorr](http://rmjarvis.github.io/TreeCorr/html/overview.html) package to do a fast weighted correlation function, using the "[optimal](http://adsabs.harvard.edu/abs/2002A%26A...386..784M)" (alpha-1) weights for magnification. 

The functions in catalog_munging.py are not too important, but just dealt with the format of the CFHTLenS catalogs that I had access to. This analysis is done using the CFHTLenS clusters (publicly available [here](https://zenodo.org/record/51291#.V5e3PZMrK9s)) as lenses, and the CFHTLenS Lyman-break galaxies (which I do not believe are yet public) as sources. This could be easily altered to apply to any background sources, as long as you have position data and an estimate of alpha(m). Since this was just a quick-and-dirty analysis, several important steps were neglected, like properly checking for consistent masking of the lenses, sources, and random catalogs.
