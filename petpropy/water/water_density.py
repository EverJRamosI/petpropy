"""
### water_density.py
#### Functions:
    - rho_water_equations
    - mccain
    
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
@email: everramosisla@gmail.com
"""

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from _tools.tools import vectorize_decorator

@vectorize_decorator
def rho_water_equation(Bw: float, *, S: int=10000) -> float:
    """This is a equation for determining the water density.

    .. math::
        rho_w = (62.4*gamma_w)/Bw
    
    Args:
        Bw (float): Water Volumetric Factor [RB/STB]
        S (int, optional): Salinity [ppm]. Defaults to 10000.

    Returns:
        rho_w (float): Water Density [lb/ft^3]
    """    
    S /= 10000
    
    gamma_w = 1 + (0.695e-6*S)
    
    rho_w = (62.4*gamma_w) / Bw
    
    return round(rho_w, 5)

@vectorize_decorator
def mccain(Bw: float, *, S: int=10000) -> float:
    """This is a correlation for determining the water density, through McCain W.D.

    Args:
        Bw (float): Water Volumetric Factor [RB/STB]
        S (int, optional): Salinity [ppm]. Defaults to 10000.

    Returns:
        rho_w (float): Water Density [lb/ft^3]
    """    
    S /= 10000
    
    rho_w1 = 62.368 + (0.438603*S) + (1.60074e-3*(S**2))
    
    rho_w = rho_w1/Bw
    
    return round(rho_w, 5)


#print('rho', rho_water_basic(1.02806, S=20000))
#print('mcain', mccain(1.02806, S=20000))
