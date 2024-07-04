"""
### gas_viscosity.py
#### Functions: 
    - lee_gonzalez_eakin

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
def lee_gonzalez_eakin(P: float|int, T: float|int, Mg: float|int, z: float) -> float:
    """This is a correlation of Lee A.L., González M.H. & Eakin B.E.\n
    For determining the gas viscosity (mu_g).

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        Mg (float | int): Weight Molecular Gas (lbs/lb-mol)
        z (float | int): Gas Compressibility Factor [dimensionless]

    #### Returns:
        mu_g (float): Gas Viscosity (cp) 
    """    
    from math import exp
    
    K = ((9.4 + 0.02 * Mg) * (T**1.5))/(209 + 19 * Mg + T)
    X = 3.5 + (986/T) + (0.01 * Mg)
    Y = 2.4 - (0.2 * X)
    
    rho_g = 1.4935e-3 * ((P * Mg)/(z * T))
    
    mu_g = round((K * exp(X*(rho_g**Y)))/(10**4), 7)
    
    return float(mu_g)


