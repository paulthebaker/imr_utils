# ringutils.py
# pythonic implementations of some BH ringdown computation tools
#  masses in Msun, freq in Hz, spins and qualities dimensionless
#
# (c) 2014 Paul T. Baker, paul.baker@ligo.org
# licence: GNU GPLv3 <http://www.gnu.org/licenses/gpl.txt>

import numpy as np

##################
# Constants used 
##################
__C = 299792458.0 # speed of light (m/s)
__GMsun = 1.32712440018e+20 # (Newton's G)*(mass of sun) (m^3/s^2)
__Tsun = __GMsun/(__C*__C*__C) # "sec per solar mass"

#######################
# fit coefficients for QNMs f
#  for l=2,3,4 and m=-l,...,l
#   From Berti et al. (2008)
#    <5% error in a=[0,0.99] 
#    note: Berti calls them 1-3 instead of 0-2
#
#  F = f0 + f1 * (1-a)^f2
#  Q = q0 + q1 * (1-a)^q2
#######################

_qnm_fs = np.array([
    [
        [None, None, None], # BLANK DATA
        [None, None, None], # BLANK DATA
        [1.5251, -1.1568, 0.1292], # l=2, m=2
        [0.6000, -0.2339, 0.4175], # l=2, m=1
        [0.4437, -0.0739, 0.3350], # l=2, m=0
        [0.3441, 0.0293, 2.0010], # l=2, m=-1
        [0.2938, 0.0782, 1.3546], # l=2, m=-2
        [None, None, None], # BLANK DATA
        [None, None, None] # BLANK DATA
    ],
    [
        [None, None, None], # BLANK DATA
        [1.8956, -1.3043, 0.1818], # l=3, m=3
        [1.1481, -0.5552, 0.3002], # l=3, m=2
        [0.8345, -0.2405, 0.4095], # l=3, m=1
        [0.6873, -0.09282, 0.3479], # l=3, m=0
        [0.5751, 0.02508, 3.1360], # l=3, m=-1
        [0.5158, 0.08195, 1.4084], # l=3, m=-2
        [0.4673, 0.1296, 1.3255], # l=3, m=-3
        [None, None, None] # BLANK DATA
    ],
    [
        [2.3000, -1.5056, 0.2244], # l=4, m=4
        [1.6869, -0.8862, 0.2822], # l=4, m=3
        [1.2702, -0.4685, 0.3835], # l=4, m=2
        [1.0507, -0.2478, 0.4348], # l=4, m=1
        [0.9175, -0.1144, 0.3511], # l=4, m=0
        [0.7908, 0.02024, 5.4628], # l=4, m=-1
        [0.7294, 0.07842, 1.5646], # l=4, m=-2
        [0.6728, 0.1338, 1.3413], # l=4, m=-3
        [0.6256, 0.1800, 1.3218] # l=4, m=-4
    ]
])


_qnm_qs = np.array([
    [
        [None, None, None], # BLANK DATA
        [None, None, None], # BLANK DATA
        [0.7000, 1.4187, -0.4990], # l=2, m=2
        [-0.3000, 2.3561, -0.2277], # l=2, m=1
        [4.0000, -1.9550, 0.1420], # l=2, m=0
        [2.0000, 0.1078, 5.0069], # l=2, m=-1
        [1.6700, 0.4192, 1.4700], # l=2, m=-2
        [None, None, None], # BLANK DATA
        [None, None, None] # BLANK DATA
    ],
    [
        [None, None, None], # BLANK DATA
        [0.9000, 2.3430, -0.4810], # l=3, m=3
        [0.8313, 2.3773, -0.3655], # l=3, m=2
        [23.8450, -20.7240, 0.03837], # l=3, m=1
        [6.7841, -3.6112, 0.09480], # l=3, m=0
        [3.0464, 0.1162, -0.2812], # l=3, m=-1
        [2.9000, 0.3356, 2.3050], # l=3, m=-2
        [2.5500, 0.6576, 1.3378], # l=3, m=-3
        [None, None, None] # BLANK DATA
    ],
    [
        [1.1929, 3.1191, -0.4825], # l=4, m=4
        [1.4812, 2.8096, -0.4271], # l=4, m=3
        [-3.6000, 7.7749, -0.1491], # l=4, m=2
        [14.0000, -9.8240, 0.09047], # l=4, m=1
        [7.0000, -2.7934, 0.1708], # l=4, m=0
        [4.6000, -0.4038, 0.4629], # l=4, m=-1
        [4.0000, 0.2777, 2.0647], # l=4, m=-2
        [3.7000, 0.5829, 1.6681], # l=4, m=-3
        [3.4000, 0.8696, 1.4074] # l=4, m=-4
    ]
], np.float64)

def f_lm(l, m):
    """get f coefficients for QNM l,m
    l=2,3,4 and -l <= m <= l
    """
    try:
        f = _qnm_fs[l-2][4-m]
        if f[0]==None:
            raise ValueError
    except IndexError:
        print("l=2,3,4; m=-l,...,+l ONLY")
    except ValueError:
        print("invalid l,m index: abs(m) > l")
    else:
        return f

def q_lm(l, m):
    """get q coefficients for QNM l,m
    l=2,3,4 and -l <= m <= l
    """
    try:
        q = _qnm_qs[l-2][4-m]
        if q[0]==None:
            raise ValueError
    except IndexError:
        print("l=2,3,4; m=-l,...,+l ONLY")
    except ValueError:
        print("invalid l,m index: abs(m) > l")
    else:
        return q

########################
# ringdown parameters:
#  F - ringdown dimensionless freq = Mw
#  Q - ringdown quality
#  f - ringdown frequency
#  M - final BH mass
#  a - final BH spin
########################

# getting F,Q from a
def FQ_from_a(a, l=2, m=2):
    """return dimensionless ringdown frequency and quality
    given spin for mode l,m; currently only l=2,3,4
    """
    f = f_lm(l, m)
    q = q_lm(l, m)
    F = f[0] + f[1] * (1.0 - a)**f[2]
    Q = q[0] + q[1] * (1.0 - a)**q[2]
    return (F, Q)

def F_from_a(a, l=2, m=2):
    """return dimensionless ringdown frequency
    given spin for mode l,m; currently only l=2,3,4
    """
    f = f_lm(l, m)
    F = f[0] + f[1] * (1.0 - a)**f[2]
    return F

def Q_from_a(a, l=2, m=2):
    """return dimensionless ringdown quality
    given spin for mode l,m; currently only l=2,3,4
    """
    q = q_lm(l, m)
    Q = q[0] + q[1] * (1.0 - a)**q[2]
    return q


# getting f,Q from M,a
def f_from_Ma(M, a, l=2, m=2):
    """return ringdown frequency given BH mass and spin
    for mode l,m; currently only l=2,3,4
    """
    F = F_from_a(a, l, m)
    f = F / (__Tsun*M) / (2.0*np.pi)
    return f

def fQ_from_Ma(M, a, l=2, m=2):
    """return ringdown freq and quality l=m=2 given BH mass and spin"""
    F, Q = FQ_from_a(a, l, m)
    f = F / (__Tsun*M) / (2.0*np.pi)
    return (f, Q)


# getting M,a from f,Q
def a_from_Q(Q, l=2, m=2):
    """return BH spin given ringdown quality for mode l,m
    """
    q = q_lm(l, m)
    a = 1.0 - np.power( (Q-q[0])/q[1], 1.0/q[2])
    return a

def Ma_from_fQ(f, Q, l=2, m=2):
    """return BH mass and spin given ringdown quality and freq
    for mode l,m; currently only l=2,3,4
    """
    a = a_from_Q(Q)
    F = F_from_a(a, l, m)
    M = F / (2.0*np.pi*f) / __Tsun
    return (M, a)

def M_from_fQ(f, Q, l=2, m=2):
    """return BH mass given ringdown quality and freq l=m=2"""
    a = a_from_Q(Q, l, m)
    F = F_from_a(a, l, m)
    M = F / (2.0*np.pi*f) / __Tsun
    return M


