"""
### gas_density.py
#### Functions:
    - rho_gas

This file is part of My Python Library.

My Python Library is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

My Python Library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with My Python Library. If not, see <http://www.gnu.org/licenses/>.

@author: Ever J. Ramos I.
"""

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from _tools.tools import vectorize_decorator

@vectorize_decorator
def rho_gas(P: float|int, T: float|int, gamma_g: float|int, z: float) -> float:
    """This is a equation of the density of gas.
    
    .. math::
        rho_g = 2.70*(P*gamma_g)/(z*T)

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_g (float | int): Specific Gravity Gas
        z (float | int): Gas Compressibility Factor [dimensionless]

    #### Returns:
        rho_g (float): Density of Gas [lbs/ft^3]
    """    
    rho_gas = round(2.70 * ((P * gamma_g)/(z * T)), 7)
    
    return float(rho_gas)



