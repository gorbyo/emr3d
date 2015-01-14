#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy.linalg import norm
from numpy import sqrt, pi


#Define global values
ng = 32                     #number of grid points 
length_system = 3.e-3       #physical length of system
nsp = 2                     #number of species
np  = 256                   #number of particles per species (np > ng)
duration = 1024             #duration of simulation
q_e_2_m_e = -1.758820150e11 #electron charge-to-mass ratio
q_e_2_m_p = 9.57883392e7    #proton charge-to-mass ratio
mass_electron = 9.10938215e-31
m_e_2_m_p = 5.4461702177e-4
m_p_2_m_e = 1836.15267247
one_amu = 1.660538782e-27
m_p = 1.672621637e-27
m_p_amu = 1.00727646677
m_4He = 6.64465620e-27
q_e = 1.602176487e-19
c_light = 299792458.
mu_0 = 4.*pi*1.e-7
eps_0 = 8.85418781762e-12
#dt = 1.

#

dt = 2.e-0/sqrt(np/length_system**3*q_e**2/eps_0/mass_electron)


def gamma_lorentz_factor(velocity):
    """return the lorentz factor"""
    speed = norm(velocity)
    if speed > c_light:
        return "speed has to be less than speed of light!"
    else:
        return 1./sqrt(1.-(speed/c_light)**2)
