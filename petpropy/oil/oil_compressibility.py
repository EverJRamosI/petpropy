"""
### oil_compressibility.py
#### Functions:
    - vazquez_beggs
    - petrosky_farshad
    - kartoatmodjo_schmidt
    - maccain_rollins_villena
    
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
def vazquez_beggs(P: float|int, T: float|int, gamma_api: float|int, gamma_gas: float, Rs: float, *, P_sp: float=0, T_sp: float=0) -> float:
    """This is a equation for determining the oil compressibility, through Vazquez M.E. & Beggs H.D.\n
    
    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]
        P_sp (float, optional): Separator Pressure [psia]. Defaults to 0.
        T_sp (float, optional): Separator Temperature [oR]. Defaults to 0.

    #### Returns:
        co (float): Oil Compressibility [psia^-1]
    """    
    from math import log10
    
    if P_sp != 0 and T_sp != 0:
        T_sp -= 460
        gamma_gc = gamma_gas * (1 + (0.1595 * (gamma_api**0.4078) * (T_sp**(-0.2466)) * log10(P_sp/114.7)))
        co = ((-1433) + (5*Rs) + (17.2*(T-460)) - (1180*gamma_gc) + (12.61*gamma_api)) / (P*(10**5))
    else:
        co = ((-1433) + (5*Rs) + (17.2*(T-460)) - (1180*gamma_gas) + (12.61*gamma_api)) / (P*(10**5))
    
    return abs(co)

@vectorize_decorator
def petrosky_farshad(P: float|int, T: float|int, gamma_api: float|int, gamma_gas: float, Rs: float) -> float:
    """This is a equation for determining the oil compressibility, through Petrosky G.E. & Farshad F.F.\n

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]

    #### Returns:
        co (float): Oil Compressibility [psia^-1]
    """    
    co = (1.705e-7) * (Rs**0.69357) * (gamma_gas**0.1885) * (gamma_api**0.3272) * ((T-460)**0.6729) * (P**(-0.5906))
    
    return abs(co)

@vectorize_decorator
def kartoatmodjo_schmidt(P: float|int, T: float|int, gamma_api: float|int, gamma_gas: float, Rs: float, *, P_sp: float=0, T_sp: float=0) -> float:
    """This is a equation for determining the oil compressibility, through Kartoatmodjo T. & Schmidt Z.\n
    
    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]
        P_sp (float, optional): Separator Pressure [psia]. Defaults to 0.
        T_sp (float, optional): Separator Temperature [oR]. Defaults to 0.

    #### Returns:
        co (float): Oil Compressibility [psia^-1]
    """  
    from math import log10
    
    if P_sp != 0 and T_sp != 0:
        gamma_gc = gamma_gas * (1 + (0.1595 * (gamma_api**0.4078) * (T_sp**(-0.2466)) * log10(P_sp/114.7)))
        co = (6.8257 * (Rs**0.5002) * (gamma_api**0.3613) * ((T-460)**0.76606) * (gamma_gc**0.35505)) / (P*(10**6))
    else:
        co = (6.8257 * (Rs**0.5002) * (gamma_api**0.3613) * ((T-460)**0.76606) * (gamma_gas**0.35505)) / (P*(10**6))
    
    return abs(co)

@vectorize_decorator
def maccain_rollins_villena(P: float|int, T: float|int, gamma_api: float|int, gamma_gas: float, *, Rs: float=0, Pb: float=0) -> float:
    """This is a equation for determining the oil compressibility, through McCain, Rollins & Villena\n

    Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        Rs (float, optional): Solution Gas-Oil Ratio [PCN/BN]. Defaults to 0.
        Pb (float, optional): Bubble Pressure [psia]. Defaults to 0.

    Returns:
        co (float): Oil Compressibility [psia^-1]
    """    
    from math import log, exp
    
    if Pb == 0 and Rs == 0:
        A = (-7.114) - (1.394*log(P)) + (0.981*log(T)) + (0.770*log(gamma_api)) + (0.446*log(gamma_gas))
        co_b = exp(A)
    elif Pb == 0:
        A = (-7.663) - (1.497*log(P)) + (1.115*log(T)) + (0.533*log(gamma_api)) + (0.184*log(Rs))
        co_b = exp(A)
    else:
        A = (-7.573) - (1.450*log(P)) - (0.383*log(Pb)) + (1.402*log(T)) + (0.256*log(gamma_api)) + (0.449*log(Rs))
        co_b = exp(A)
        
    return abs(co_b)



#print(c_oil(4000, (180+460), 31, 0.95, Rs=760, model='VB'))
#print('vazquez', vazquez_beggs(4000, (180+460), 582, 31, 0.95))
#print('petrosky', petrosky_farshad(4000, (180+460), 624, 31, 0.95))
#print('kartoatmo', kartoatmodjo_schmidt(4000, (180+460), 555, 31, 0.95))
#print('maccain', maccain_rollins_villena(2000, (180+460), 31, 0.95, Rs=516, Pb=2500))