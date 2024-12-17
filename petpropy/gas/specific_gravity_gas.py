"""
### specific_gravity_gas.py
#### Functions:
    - gamma_gas
    
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

def gamma_gas(m_gas: float) -> float:
    """This is a equation for determining the Specific Gravity Gas.
    
    .. math::
        gamma_g = m_g/28.96

    #### Args:
        m_gas (float): Weight Molecular Gas [lb/lb-mol]

    #### Returns:
        gamma_g (float): Specific Gravity Gas [fraction]
    """    
    gamma_gas = m_gas / 28.96
    return gamma_gas

