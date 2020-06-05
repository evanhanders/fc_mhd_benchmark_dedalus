import numpy as np

import dedalus.public as de

class CartesianBoxCF:
    """
    Creates a dedalus domain which is a cartesian box.

    The bases are Fourier in the horizontal directions and Chebyshev in the vertical.

    # Attributes
        bases (list) : A list of Dedalus bases
        domain (Domain) : The Dedalus domain
    """

    def __init__(self, resolution=(256, 64), basis_names=('x', 'z'), height=1, aspect=4, dealias=3/2, grid_dtype=np.float64, mesh=None):
        """

        # Arguments
            resolution (tuple, optional) :
                A tuple of ints of the basis coefficient resolutions (horizontal, then vertical)
            basis_names (tuple, optional) :
                A tuple of strings containing the desired names of the bases (horizontal, then vertical)
            height (float, optional) :
                The height of the atmosphere, which runs from [0, height]
            aspect (float, optional) :
                The aspect ratio of the atmosphere. Horizontally, the box runs from [-height*aspect/2, height*aspect/2]
            dealias (float, optional) :
                The spectral dealiasing factor
            grid_dtype (datatype, optional) :
                The datatype of Dedalus in grid space
            mesh (tuple, optional) :
                The distribution mesh of the processors for running in parallel
        """
        self.bases = []
        for name, resolution in zip(basis_names, resolution):
            if name == basis_names[-1]:
                self.bases.append(de.Chebyshev(name, resolution, interval=[0, height], dealias=dealias))
            else:
                self.bases.append(de.Fourier(name, resolution, interval=[-height*aspect/2, height*aspect/2], dealias=dealias))
        self.domain = de.Domain(self.bases, grid_dtype=grid_dtype, mesh=mesh)
