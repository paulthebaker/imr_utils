# inspiral.py
# pythonic implementations of some BH inspiral computation tools
#  masses in Msun, freq in Hz, spins and qualities dimensionless
#
# (c) 2014 Paul T. Baker, bakerp@geneseo.edu
# licence: GNU GPLv3 <http://www.gnu.org/licenses/gpl.txt>

import numpy as np

##############################
# binary orbital parameters:
#  m1 - component mass 1
#  m2 - component mass 2
#  M - system total mass
#  q - mass ratio (>1)
#  eta - symmetric mass ratio
##############################

def Mq_from_m1m2(m1, m2):
    """return total mass and ratio given component masses"""
    M = m1 + m2
    q = max(m1/m2, m2/m1)
    return (M, q)

def m1m2_from_Mq(M, q):
    """return component masses given total mass and ratio
    invert:  M=m1+m2 and q=m2/m1
    """
    m1 = M/(q+1.0)
    m2 = q*m1
    return (m1, m2)


def eta_from_m1m2(m1, m2):
    """return symmetric mass ratio given component masses"""
    return m1*m2 / (m1+m2)**2


def eta_from_Mq(M, q):
    """return symmetric mass ratio given total mass and ratio"""
    m1, m2 = m1m2_from_Mq(M, q)
    return m1*m2 / (m1+m2)**2


###########################
# final BH parameters: for initially non-spinning BBH
#  Mf - final BH mass... different than Mtot due to radiation
#  a - final BH spin
#  eta - initial symmetric mass ratio
#  m1 - initial component mass 1
#  m2 - initial component mass 2
#  Mi - initial system total mass
#  q - initial mass ratio
#  e - fractional energy radiated
###########################

def a_from_eta(eta): 
    """return BH spin given symmetric mass ratio
    only good for initially non-spinning BBH
    """
    return np.sqrt(12.0)*eta - 2.9*eta*eta


def a_from_m1m2(m1, m2): 
    """return BH spin given component masses
    only good for initially non-spinning BBH
    """
    eta = eta_from_m1m2(m1, m2)
    return a_from_eta(eta)


def a_from_Miq(Mi, q): 
    """return BH spin given initial Mtot and mass ratio
    only good for initially non-spinning BBH
    """
    eta = eta_from_Mq(Mi, q)
    return a_from_eta(eta)


def Mf_from_etaMi(eta, Mi):
    """return final BH mass given initial Mtot and symmetric mass ratio
    only good for initially non-spinning BBH
    """
    return ( 1 + (np.sqrt(8.0/9.0) - 1)*eta - 0.498*eta*eta) * Mi


def Mf_from_m1m2(m1, m2): 
    """return final BH mass given component masses
    only good for initially non-spinning BBH
    """
    eta = eta_from_m1m2(m1, m2)
    Mi = m1 + m2
    return Mf_from_etaMi(eta, Mi)


def Mf_from_Miq(Mi, q): 
    """return final BH mass given initial Mtot and mass ratio
    only good for initially non-spinning BBH
    """
    eta = eta_from_Mq(Mi, q)
    return Mf_from_etaMi(eta, Mi)


