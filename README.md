# imr_utils

Python utilities for computing binary black hole merger paremters including QNMs


## sub modules

inspiral.py
 * convert between component masses, total mass and mass ratio (or symmetric mass ratio)
 * compute final black hole mass and spin from any of the above

ringdown.py
 * compute black hole ringdown QNM frequency and quality factor
 * handles all n=1, l=2,3,4 modes (m=[-l,+l])
 * compute real or dimensionless frequency from mass and spin
 * compute BH mass and spin from frequency and quality
 

QNM calculations based on approximation of [Berti et al. (2008)](http://arxiv.org/abs/gr-qc/0512160)

