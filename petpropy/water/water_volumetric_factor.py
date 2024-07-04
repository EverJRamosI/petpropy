"""
### water_volumetric_factor.py
#### Functions:
    - mccain
    - mccoy
    - factor_bw
    
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
def mccain(P: float|int, T: float|int) -> float:
    """This is a correlation for determining the water volumetric factor, through McCain W.D.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]

    #### Returns:
        Bw (float): Water Volumetric Factor [RB/STB]
    """    
    dV_wT = (-1.0001e-2) + (1.33391e-4*(T-460)) + (5.50654e-7*((T-460)**2))
    dV_wP = (-1.95301e-9*P*(T-460)) - (1.72834e-13*(P**2)*(T-460)) - (3.58922e-7*P) - (2.25341e-10*(P**2))
    
    Bw = (1 + dV_wP)*(1 + dV_wT)
    
    return Bw

@vectorize_decorator
def mccoy(P: float|int, T: float|int, *, S: int=10000, type_gas: int=1) -> float:
    """This is a correlation for determining the water volumetric factor, through McCoy R.L.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        S (int, optional): Salinity [ppm]. Defaults to 10000.
        type_gas (int, optional): Type Gas: gas free=1 | gas saturated=2. Defaults to 1.

    #### Returns:
        Bw (float): Water Volumetric Factor [RB/STB]
    """    
    
    T -= 460
    S /= 10000
        
    if type_gas == 1:
        
        A = 0.9947 + (5.8e-6*T) + (1.02e-6*(T**2))
        B = (-4.228e-6) + (1.8376e-8*T) - (6.77e-11*(T**2))
        C = (1.3e-10) - (1.3855e-12*T) + (4.285e-15*(T**2))
        
        B_wp = A + (B*P) + (C*(P**2))
        
    elif type_gas == 2:
        
        A = 0.9911 + (6.35e-5*T) + (8.5e-7*(T**2))
        B = (-1.093e-6) - (3.497e-9*T) + (4.57e-12*(T**2))
        C = (-5.0e-11) + (6.429e-13*T) - (1.43e-15*(T**2))
        
        B_wp = A + (B*P) + (C*(P**2))
        
    Bw = (1 + (S * ((5.1e-8*P) + ((5.47e-6-(1.95e-10*P))*(T-60))-((3.23e-8-(8.5e-13*P))*((T-60)**2))))) * B_wp
        
    return Bw

@vectorize_decorator
def factor_bw(P: float|int, Pb: float|int, Bwb: float, cw: float) -> float:
    """This is a equation for determining the water volumetric factor.
    
    .. math::
        Bw = Bwb*exp(cw*(Pb-P))

    #### Args:
        P (float | int): Pressure [psia]
        Pb (float | int): Bubble Pressure [psia]
        Bwb (float): Water Volumetric Factor at the Bubble Point [RB/STB]
        cw (float): Water Compressibility [psia^-1]

    #### Returns:
        Bw (float): Water Volumetric Factor [RB/STB]
    """    
    from math import exp
    
    Bw = Bwb * exp(cw*(Pb-P))
    
    return Bw


#print('mccain', mccain(5000, (200+460)))
#print('mccoy', mccoy(5000, (200+460), S=20000, type_gas=2))

# import matplotlib.pyplot as plt

# presiones = list(range(500, 5100, 100))

# temperatura = 200 + 460

# #factor = [mccain(p, temperatura) for p in presiones]
# factor = [mccoy(p, temperatura, S=20000, type_gas=2) for p in presiones]

# print(factor)

# plt.plot(presiones, factor)
# plt.show()