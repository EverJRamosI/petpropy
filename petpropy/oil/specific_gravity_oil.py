"""
### specific_gravity_oil.py
#### Functions:
    - gamma_oil
    - api
    
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

def gamma_oil(rho_oil: float, rho_water: float) -> float:
    """This is a equation for determining the specific Gravity Oil.
    
    .. math::
        gamma_o = rho_o/rho_w

    Args:
        rho_oil (float): Oil Density [lb/ft^3]
        rho_water (float): Water Density [lb/ft^3]

    Returns:
        gamma_o (float): Specific Gravity Oil
    """    
    gamma_o = rho_oil / rho_water
    return gamma_o

def api(gamma_oil: float) -> float:
    """This is equation for determining the API gravity.
    
    .. math::
        API = (141.5/gamma_o)-131.5

    Args:
        gamma_oil (float): Specific Gravity Oil

    Returns:
        API (float): API Gravity
    """    
    gamma_api = (141.5 / gamma_oil) - 131.5
    return gamma_api

