import numpy as np

from fce_dedalus.atmospheres.atmosphere_base import IdealGasAtmosphere


class Polytrope(IdealGasAtmosphere):
    """
    An atmosphere with a simple polytropic stratification.

    Given a polytropic index, m, and a number of density scale heights, nρ, the stratification is:

        T  = 1 + (∇T)z
        ∇T = (e^{-nρ/m} - 1)/Lz, a constant
        ρ  = T^m

    with
        g = -(1 + m)∇T

    and
        Lz = 1.
    """

    def __init__(self, nρ=3, ε=0.5, ɣ=5./3, *args, **kwargs):
        """
        Initialize the polytrope

        # Arguments
            nρ (float) : The number of density scale heights
            ε (float)  : The superadiabatic excess
            ɣ (float)  : The adiabatic index
            args, kwargs :
                Additional arguments and keyword arguments to pass through to the parent class init function.
        """
        self['nρ']   = nρ
        self['ε']    = ε
        self['m_ad'] = 1/(ɣ - 1)
        self['m']    = self['m_ad'] - self['ε']
        super(Polytrope, self).__init__(*args, ɣ=ɣ, **kwargs)

    def _set_stratification(self):
        """ Creates functions for the atmospheric density and temperature stratification. """
        self['L']  = 1
        self['∇T'] = (np.exp(-self['nρ']/self['m']) - 1)/self['L']
        self['T']  = lambda z : 1 + self['∇T']*z
        self['ρ']  = lambda z : self['T'](z)**self['m']
