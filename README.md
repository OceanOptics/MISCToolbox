MISC Lab Toolbox
================

**Matlab tools for oceanographic analysis focusing on the inherent optical properties (IOPs) of the open ocean.**

Features
--------
The functions available are:
* compute_bbp.m: Compute the particulate backscattering b_bp from the backscattering bb
* correct_npq.m: Correct for non photochemical quenching using Xing et al. (2012) and/or Sackmann et al. (2008)
* etimate_mld.m: Estimate the mixed layer depth (MLD) with one of the following method: fixed temperature threshold, fixed density threshold, variable density threshold or fixed density gradient.
* meshprofile.m: Interpolate data between profile (often used with scatter3m)
* need_npqc.m: Determine if need a non photochemical quenching
* scatter3m.m: 4D visualization with earth map (latitude, longitude, depth and measure)

Requirements
------------
To work properly the toolbox need those features in matlab path (addpath):
* betasw_ZHH2009.m from Xiaodong Zhang available [here](https://github.com/ooici/ion-functions/blob/master/ion_functions/data/matlab_scripts/flort/betasw_ZHH2009.m)
* gsw_matlab_v3_04 from TEOS-10 available [here](https://github.com/TEOS-10/GSW-Matlab/releases)
* lr2.m a robust linear regression type II, can substitute it by regress() instead
