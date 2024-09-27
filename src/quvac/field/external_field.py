'''
This script provides uniform ExternalField class to unite all
participating fields in one interface
'''
import os

import numpy as np

from quvac.field.abc import Field
from quvac.field.paraxial_gaussian import ParaxialGaussianAnalytic
from quvac.field.maxwell import ParaxialGaussianMaxwell


class ExternalField(Field):
    '''
    Class to unite several participating fields under
    one interface

    Parameters
    ----------
    fields_params: list of dicts (field_params)
        External fields
    grid: (1d-np.array, 1d-np.array, 1d-np.array)
        xyz spatial grid to calculate fields on 
    '''
    def __init__(self, fields_params, grid, nthreads=None):
        self.fields = []
        self.grid = grid

        self.nthreads = nthreads if nthreads else os.cpu_count()

        for field_params in fields_params:
            self.setup_field(field_params)

    def setup_field(self, field_params):
        field_type = field_params["field_type"]
        match field_type:
            case "paraxial_gaussian_analytic":
                field = ParaxialGaussianAnalytic(field_params, self.grid)
            case "paraxial_gaussian_maxwell":
                field = ParaxialGaussianMaxwell(field_params, self.grid, nthreads=self.nthreads)
            case _:
                raise NotImplementedError(f"We do not support '{field_type}' field type")
        self.fields.append(field)
            
    def calculate_field(self, t, E_out=None, B_out=None):
        E_out_= [np.zeros(E_out[0].shape, dtype=np.complex128) for _ in range(3)]
        B_out_= [np.zeros(E_out[0].shape, dtype=np.complex128) for _ in range(3)]
        for field in self.fields:
            field.calculate_field(t, E_out=E_out_, B_out=B_out_)
        for i in range(3):
            E_out[i] += np.real(E_out_[i])
            B_out[i] += np.real(B_out_[i])

