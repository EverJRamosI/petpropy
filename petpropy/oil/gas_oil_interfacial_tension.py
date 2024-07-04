"""
### gas_oil_interfacial_tension.py
#### Functions:
    - baker_swedloff
    
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
def baker_swedloff(P: float|int, T: float|int, gamma_api: float|int) -> float:
    """Correlation of Gas-Oil Interfacial Tension, through Baker O. & Swerdloff W.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity

    #### Returns:
        sigma_o (float): Gas-Oil Interfacial Tension [dynes/cm]
    """    
    F = 1 - (0.024*(P**0.45))
    
    if (T-460) <= 68:
        sigma_68 = 39 - (0.2571*gamma_api)
        return sigma_68*F
    elif (T-460) >= 100: 
        sigma_100 = 37.5 - (0.2571*gamma_api)
        return sigma_100*F
    else:
        sigma_68 = 39 - (0.2571*gamma_api)
        sigma_100 = 37.5 - (0.2571*gamma_api)
        sigma_T = sigma_68 - ((((T-460) - 68)*(sigma_68-sigma_100))/32)
        return sigma_T*F
    
#import matplotlib.pyplot as plt

#presiones = list(range(500, 3000, 100))

# temperatura = 90+460

# api_oil = 31

# tension = [baker_swedloff(p, temperatura, api_oil) for p in presiones]

#print(presiones)
#print(tension)

#plt.plot(tension, presiones)
#plt.show()

