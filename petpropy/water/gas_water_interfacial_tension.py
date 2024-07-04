"""
### gas_water_interfacial_tension.py
#### Functions
    - jennings_newman
    
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
def jennings_newman(P: float|int, T: float|int) -> float:
    """This is a correlation for determining the gas-water interfacial tension, through Jennings H.Y. & Newman G.H.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]

    #### Returns:
        sigma_gw (float): Gas-Water Interfacial Tension [dynes/cm]
    """    
    T -= 460
    
    A = 79.1618 - (0.118978*T)
    B = (-5.28473e-3) + (9.87913e-6*T)
    C = (2.33814 - (4.57194e-4*T) - (7.52678e-6*(T**2))) * (10**(-7))
    
    sigma_gw = A + (B*P) + (C*(P**2))
    
    return round(sigma_gw, 7)

#ps = list(range(500, 3000, 500))
#it = jennings_newman(ps, (200+460))
#print('jenning', it)