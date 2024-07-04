"""
### water_viscosity.py
#### Functions:
    - van_wingen
    - matthews_russel
    - mccain
    - mccoy
    
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
def van_wingen(T: float|int) -> float:
    """This is a correlation for determining the water viscosity, through Van Wingen N.

    #### Args:
        T (float | int): Temperature [oR]

    #### Returns:
        mu_w (float): Water Viscosity [cp]
    """    
    from math import exp
    
    T -= 460
    
    mu_w = exp(1.003 - (1.479e-2*T) + (1.982e-5*(T**2)))
    
    return round(mu_w, 6)

@vectorize_decorator
def matthews_russel(P: float|int, T: float|int, *, S: int=10000) -> float:
    """This is a correlation for determining the water viscosity, through Matthews C.S. & Russel D.G.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        S (int, optional): Salinity [ppm]. Defaults to 10000.

    #### Returns:
        mu_w (float): Water Viscosity [cp]
    """    
    T -= 460
    S /= 10000
    
    A = (-0.04518) + (0.009313*S) - (0.000393*(S**2))
    B = 70.634 + (0.09576*(S**2))
    
    mu_w1 = A + (B/T)
    f = 1 + (3.5e-12*(P**2)*(T-40))
    
    mu_w = mu_w1*f
    
    return round(mu_w, 6)

@vectorize_decorator
def mccain(P: float|int, T: float|int, *, S: int=10000) -> float:
    """This is a correlation for determining the water viscosity, through McCain W.D.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        S (int, optional): Salinity [ppm]. Defaults to 10000.

    #### Returns:
        mu_w (float): Water Viscosity [cp]
    """
    T -= 460
    S /= 10000
    
    A = 109.574 - (8.40564*S) + (0.313314*(S**2)) + (8.72213e-3*(S**3))
    B = (-1.12166) + (2.63951e-2*S) - (6.79461e-4*(S**2)) - (5.47119e-5*(S**3)) + (1.55586e-6*(S**4))
    
    mu_w1 = A*(T**B)
    
    mu_w = (0.9994 + (4.0295e-5*P) + (3.1062e-9*(P**2))) * mu_w1
    
    return round(mu_w, 6)

@vectorize_decorator
def mccoy(T: float|int, *, S: int=10000) -> float:
    """This is a correlation for determining the water viscosity, through McCoy R.L.

    #### Args:
        T (float | int): Temperature [oR]
        S (int, optional): Salinity [ppm]. Defaults to 10000.

    #### Returns:
        mu_w (float): Water Viscosity [cp]
    """    
    T -= 460
    S /= 10000
    T = ((5/9)*T) + 255.37
    
    mu_wp = 0.02414 * (10**(247.8/(T-140)))
    
    mu_w = (1 - (1.87e-3*(S**0.5)) + (2.18e-4*(S**2.5)) + ((T**0.5)-(1.35e-2*T))*((2.76e-3*S)-(3.44e-4*(S**1.5)))) * mu_wp
    
    return round(mu_w, 6)

# temperature = [150, 250, 350, 400]
# t = [t+460 for t in temperature]
# print('van', van_wingen((t)))
# p = [2500, 3000, 3500, 4000]
# print('mathe', matthews_russel(p, t, S=20000))
#print('van', van_wingen((200+460)))
#print('mathe', matthews_russel(5000, (200+460), S=20000))
#print('mcain', mccain(5000, (200+460), S=20000))
#print('mcoy', mccoy((200+460), S=20000))


# import matplotlib.pyplot as plt

# presiones = list(range(500, 5100, 100))

# temperatura = 200 + 460

# #viscosidad = [matthews_russel(p, temperatura, S=20000) for p in presiones]
# viscosidad = [mccain(p, temperatura, S=20000) for p in presiones]

# print(presiones)
# print(viscosidad)

# plt.plot(presiones, viscosidad)
# plt.show()