# 2D Model Waleffe Flow (MWF)
This code solves a two-vertical-mode truncation of the [Waleffe Flow (WF)]{https://doi.org/10.1017/jfm.2016.92} system and its turbulence closures, resulting in a six variable, two-dimensional system of partial differential equations, which are solved using [Dedalus]{https://dedalus-project.org/} (currently just [v2_master]{https://github.com/DedalusProject/dedalus/tree/v2_master}). See Benavides, Tuckerman, & Barkley (arxiv link coming soon...) for details on the truncation, closure choices, and use of the model.

This repository contains two directories:

1. A sample run directory, `Lx224Lz100_Re75_tube_minimal`, which contains a parameter file `params.py`, and the main code file `main.py`, as well as a sample jobscript fle, Dedalus configuration file, and a script to clean the directory used in restarting runs from zero.

2. A post-processing directory, `postproc`, which contains a movie making script and a Jupyter Lab notebook with scripts for plotting time series and snapshots from the simulations.

If you have any questions regarding the code or its implementation, don't hesitate to reach out to Santiago.Benavides@warwick.ac.uk 
