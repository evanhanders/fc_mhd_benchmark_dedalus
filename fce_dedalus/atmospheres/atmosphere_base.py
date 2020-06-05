from collections import OrderedDict

import numpy as np

def NCC(domain, field_function, constant_bases=['x',], scales=1, return_scales=1):
    """
    Turns a given function into a Dedalus non-constant coefficient (NCC).

    # Arguments
        domain (Dedalus domain) :
            The Dedalus domain object on which the problem is being solved.
        field_function (Numpy array) :
            A 1-D numpy array containing the NCC data
        constant_bases (list, optional) :
            The names of the bases over which the NCC does NOT vary (horizontal bases)
        scales (float, optional) :
            The number of data points in field_function divided by the spectral resolution of the problem.
            To ensure an NCC only expands to, e.g., 4 modes, make sure field_function only has 4 data points and scales=4/nz.
        return_scales (float, optional) :
            The scale at which to return the NCC
    """
    NCC = domain.new_field()
    NCC.meta[tuple(constant_bases)]['constant'] = True
    NCC.set_scales(scales, keep_data=False)
    NCC['g'] = field_function
    NCC.set_scales(return_scales, keep_data=False)
    return NCC


class IdealGasAtmosphere():
    """
    An abstract class defining the basic functionality of an ideal gas atmosphere for convection simulations.

    # Public Methods
    - __init__()
    - build_atmosphere()

    # Attributes
        self.atmosphere_params (OrderedDictionary) :
            A dictionary that defines the fundamentals elements of the atmosphere (T, ρ, etc)
    """
    atmosphere_params = OrderedDict()
    function_fields   = ['T', 'ρ', 'P']

    def __init__(self, *args, R=1, ɣ=5./3, **kwargs):
        """
        Initializes the class.

        # Arguments
            R (float) :
                The ideal gas constant, P = R ρ T
            ɣ (float) :
                The adiabatic index of the ideal gas (ɣ = 5/3 for monatomic, 7/5 for diatomic)
            args, kwargs :
                Additional arguments and keyword arguments for the _set_stratification() function
        """
        # Basic thermodynamics
        self['R'] = R
        self['ɣ'] = ɣ
        self['Cp'] = self['R'] * self['ɣ'] / (self['ɣ'] - 1)
        self['Cv'] = self['Cp'] = self['R']
       
        self._set_stratification(*args, **kwargs)
        self['P'] = lambda z : self['R']*( self['ρ'](z) * self['T'](z) )

    def __getitem__(self, key):
        return self.atmosphere_params[key]

    def __setitem__(self, key, value):
        self.atmosphere_params[key] = value

    def _set_stratification(self):
        """ Abstract function; should set T, ρ, and L in child classes. """
        pass

    def build_atmosphere(self, domain, problem, delta_s_bounds=None, **kwargs):
        """
        Build an atmosphere and load it into dedalus. 

        Builds non-constant-coefficient fields for T, ρ, lnρ, dz(T), dz(lnρ), and s/Cp, according to the
        prescription specified by the class in the self._set_stratification() function.

        # Arguments
            domain (Dedalus domain) :
                The Dedalus domain object for the problem being solved
            problem (Dedalus problem) :
                The Dedalus problem object for the IVP, EVP, BVP, etc. being solved.
            delta_s_bounds (list, optional) :
                A 2-element list containing the z-coordinates of the bottom and top of the region over which to calculate Δs (e.g., the convection zone)
            kwargs (dict) :
                Optional keyword arguments for the NCC() function.
        """
        z_basis = domain.bases[-1]
        z_name  = domain.bases[-1].name
        z_de    = domain.grid(-1, scales=domain.dealias)

        self['T_field']   = NCC(domain, self['T'](z_de),         scales=domain.dealias, return_scales=domain.dealias, **kwargs)
        self['ρ_field']   = NCC(domain, self['ρ'](z_de),         scales=domain.dealias, return_scales=domain.dealias, **kwargs)
        self['lnρ_field'] = NCC(domain, np.log(self['ρ'](z_de)), scales=domain.dealias, return_scales=domain.dealias, **kwargs)
        
        self['Tz_field']   = NCC(domain, 0, **kwargs)
        self['lnρz_field'] = NCC(domain, 0, **kwargs)
        self['T_field'].differentiate(z_name,   out=self['Tz_field'])
        self['lnρ_field'].differentiate(z_name, out=self['lnρz_field'])

        self['s_div_Cp'] = (1/self['ɣ']) * (self['T_field']['g'] - (self['ɣ']-1)*self['lnρ_field']['g'] )
        self['s_div_Cp'] = NCC(domain, self['s_div_Cp'], scales=domain.dealias, return_scales=domain.dealias, **kwargs)
        if delta_s_bounds is None:
            self['Δs/Cp'] = np.mean(self['s_div_Cp'].integrate(z_name)['g'])
        else:
            self['Δs/Cp'] = np.mean(basis.Interpolate(self['s_div_Cp'], delta_s_bounds[-1])['g']) \
                           -np.mean(basis.Interpolate(self['s_div_Cp'], delta_s_bounds[0])['g'])

        for k in self.atmosphere_params.keys():
            if k in self.function_fields:
                continue
            problem.parameters[k] = self[k]
            print(type(self[k]), type(self[k]) == float, type(self[k]) == int, k)
