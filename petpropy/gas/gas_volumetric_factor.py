"""
### gas_volumetric_factor.py
#### Functions:
    - B_g
    - E_g

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
def B_g(P: float|int, T: float|int, z: float, *, units: bool=False) -> float:
    """Equation of gas volumetric factor (Bg).

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        z (float): Gas Compressibility Factor [dimensionless]
        units (bool, optional): Units to select. Defaults to False for [PCY/PCN] and True for [BY/PCN].

    #### Returns:
        Bg (float): Gas Volumetric Factor [PCY/PCN | BY/PCN].
    """    
    if units:
        Bg = 0.00503 * ((z*T) / P)
    else:
        Bg = 0.02827 * ((z*T) / P)
    
    return Bg

@vectorize_decorator
def E_g(P: float|int=None, T: float|int=None, z: float=None, *, Bg: float=None, units: bool=False) -> float:
    """Equation of gas expansion factor (Eg).
    If you only want to input Bg is a kwargs.
    
    #### Args:
        P (float | int, optional): Pressure [psia]. Defaults to None.
        T (float | int, optional): Temperature [oR]. Defaults to None.
        z (float, optional): Gas Compressibility Factor [dimensionless]. Defaults to None.
        Bg (float, optional): Gas Volumetric Factor [PCY/PCN | BY/PCN]. Defaults to None.
        units (bool, optional): Units to select. Defaults to False for [PCY/PCN] and True for [BY/PCN].

    #### Returns:
        Eg (float): Gas Expansion Factor [PCN/PCY | PCN/BY].
    """    
    if Bg is not None:
        Eg = 1/Bg
    elif (P, T, z) is not None:
        if units:
            Eg = 198.8 * (P/(z*T))
        else:
            Eg = 35.37 * (P/(z*T))
    
    return Eg


# import matplotlib.pyplot as plt
# from gas_compressibility_factor import z_gas

# presiones = list(range(500, 7500, 500))
# temperatura = 654
# prpc = 680
# tepc = 485.9

# zf = [z_gas(p, temperatura, prpc, tepc) for p in presiones]
# bg = [B_g(p, temperatura, z) for p, z in zip(presiones, zf)]


# print(presiones)
# print(zf)
# print(bg)

# plt.plot(presiones, bg)
# plt.show()