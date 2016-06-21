# imr_utils

Python utilities for computing binary black hole merger paremters including QNMs

Unit conventions:
 * mass (_M_,_m_) - solar mass
 * spin (_a_) - dimensionless (ang. mom. per _M_<sup>2</sup>)
 * real frequency (_f_) - hz
 * dimensionless frequency (_F_) - dimensionless...
 * quality factor (_Q_)- dimensionless

## sub modules

inspiral.py
 * convert between component masses, total mass and mass ratio (or symmetric mass ratio)
 * compute final black hole mass and spin from any of the above

ringdown.py
 * compute black hole ringdown QNM frequency and quality factor
 * handles all _l_=2,3,4 modes for all _m_=[-_l_,+_l_], _n_=1 only
 * compute real (hz) or dimensionless frequency from mass and spin
 * compute BH mass and spin from frequency and quality
 

QNM calculations based on approximation of [Berti et al. (2005)](http://arxiv.org/abs/gr-qc/0512160)

The _l_=_m_=2 functionality is implemented in [pycbc](https://github.com/ligo-cbc/pycbc), the other modes are not.
