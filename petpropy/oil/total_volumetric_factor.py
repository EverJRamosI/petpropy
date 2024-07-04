"""
### total_volumetric_factor.py
#### Functions:
    - factor_total
    - glaso
    - almarhoun
    
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
def factor_total(Bo: float, Bg: float, Rsi: float, Rs: float) -> float:
    """This is a equation for determining of the factor volumetric total.
    
    .. math::
        Bt = Bo*((Rsi-Rs)*Bg)

    Args:
        Bo (float): Oil Volumetric Factor [RB/STB]
        Bg (float): Gas Volumetric Factor [RB/PCN]
        Rsi (float): Solution Gas-Oil Ratio Initial [PCN/STB]
        Rs (float): Solution Gas-Oil Ratio [PCN/STB]

    Returns:
        Bt (float): Oil Volumetric Factor Total [RB/STB]
    """    
    Bt = Bo + ((Rsi-Rs)*Bg)
    
    return Bt

@vectorize_decorator
def glaso(P: float|int, T: float|int, Rs: float, gamma_api: float|int, gamma_gas: float) -> float:
    """This is correlation for determining the factor volumetric total, through Glaso O.

    Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        Rs (float): Solution Gas-Oil Ratio [PCN/STB]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas

    Returns:
        Bt (float): Oil Volumetric Factor Total [RB/STB]
    """    
    from math import log10
    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    F = Rs * (((T-460)**0.5)/(gamma_gas**0.3)) * (P**(-1.1089)) * (gamma_oil**(2.9*(10**((-0.00027)*Rs))))
    Bt = 10**((8.0135e-2) + (4.7257e-1*log10(F)) + (1.7351e-1*((log10(F))**2)))
    
    return Bt

def almarhoun(P: float|int, T: float|int, Rs: float, gamma_api: float|int, gamma_gas: float) -> float:
    """This is correlation for determining the factor volumetric total, through Al-Marhoun M.A.

    Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        Rs (float): Solution Gas-Oil Ratio [PCN/STB]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas

    Returns:
        Bt (float): Oil Volumetric Factor Total [RB/STB]
    """    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    F = (Rs**0.644516) * (gamma_gas**(-1.07934)) * (gamma_oil**0.724874) * (P**(-0.76191)) * (T**2.00621)
    Bt = 0.314693 + (0.106253e-4*F) + (0.18883e-10*(F**2))
    
    return Bt

#print("glaso", glaso_correlation(2000, (180+460), 433, 31, 0.95))
#print("almarhoun", almarhoun_correlation(2000, (180+460), 615, 31, 0.95))


# import matplotlib.pyplot as plt
# from dissolved_gas_oil_ratio import standing as rsst, almarhoun as rsalm

# presiones = list(range(500, 3100, 100))
# temperature = 180 + 460
# api_oil = 31
# g_gas = 0.95
# burbuja = 2500

# #disol = [rsst(p, temperature, g_gas, api_oil, Pb=burbuja) for p in presiones]
# disol = [rsalm(p, temperature, g_gas, api_oil, Pb=burbuja) for p in presiones]
# #bt = [glaso_correlation(p, temperature, rs, api_oil, g_gas) for p, rs in zip(presiones, disol)]
# bt = [almarhoun_correlation(p, temperature, rs, api_oil, g_gas) for p, rs in zip(presiones, disol)]

# print(presiones)
# print(disol)
# print(bt)

# plt.figure(figsize=(10, 6))
# plt.plot(presiones, bt, marker='o', linestyle='--')

# plt.show()
# plt.grid(True)
