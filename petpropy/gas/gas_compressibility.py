"""
### gas_compressibility.py

#### Functions:
    - papay
    - brill_beggs
    - gopal

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
def papay(P: float, T: float, Ppc: float, Tpc: float, z: float) -> float:
    """This is a correlation for determining the gas compressibility (c_gas) for method of Papay J.

    #### Args:
        P (float): Pressure [psia]
        T (float): Temperature [oR]
        Ppc (float): Pressure pseudocritical [psia]
        Tpc (float): Temperature pseudocritica [oR]
        z (float): Gas Compressibility Factor [dimensionless]

    #### Returns:
        c_gas (float): Gas Compressibility [psia^-1]
    """    
    P_pr = P / Ppc
    T_pr = T / Tpc
    
    dZr = - ((3.52) / (10**(0.9813*T_pr))) + ((0.548*P_pr) / (10**(0.8157*T_pr)))
    
    cr = (1/P_pr) - ((1/z)*dZr)
    
    cg = cr/Ppc
    
    return cg

@vectorize_decorator
def brill_beggs(P: float, T: float, Ppc: float, Tpc: float, z: float) -> float:
    """This is a correlation for determining the gas compressibility (c_gas) for method of Brill J.P. & Beggs H.D.

    #### Args:
        P (float): Pressure [psia]
        T (float): Temperature [oR]
        Ppc (float): Pressure pseudocritical [psia]
        Tpc (float): Temperature pseudocritica [oR]
        z (float): Gas Compressibility Factor [dimensionless]

    #### Returns:
        c_gas (float): Gas Compressibility [psia^-1]
    """    
    from math import exp, log
    
    P_pr = P / Ppc
    T_pr = T / Tpc
    
    A = 1.39 * (T_pr-0.92)**0.5 - 0.36 * T_pr - 0.10
    B = ((0.62-(0.23*T_pr))*P_pr) + (((0.066/(T_pr-0.86))-0.037)*(P_pr**2)) + ((0.32/(10**(9*(T_pr-1))))*(P_pr**6))
    C = 0.132 - (0.32 * log(T_pr))
    D = 10**(0.3106-(0.49*T_pr)+0.1824*(T_pr**2))
    
    dZr = ((1-A) / (((0.62-0.23*T_pr) + ((((0.132)/(T_pr-0.86))-0.074)*P_pr) + (((1.92)/(10**(9*(T_pr-1))))*(P_pr**5)))*(exp(B)))) + (C*D*(P_pr**(D-1)))
    
    cr = (1/P_pr) - ((1/z)*dZr)
    
    cg = cr/Ppc
    
    return cg

@vectorize_decorator
def gopal(P: float, T: float, Ppc: float, Tpc: float, z: float) -> float:
    """This is a correlation for determining the gas compressibility (c_gas) for method of Gopal V.N.

    #### Args:
        P (float): Pressure [psia]
        T (float): Temperature [oR]
        Ppc (float): Pressure pseudocritical [psia]
        Tpc (float): Temperature pseudocritica [oR]
        z (float): Gas Compressibility Factor [dimensionless]

    #### Returns:
        c_gas (float): Gas Compressibility [psia^-1]
    """
    P_pr = P / Ppc
    T_pr = T / Tpc
    
    if 0.2 <= P_pr <= 1.2:
        if 1.05 <= T_pr <= 1.2:
            dZr = (1.6643 * T_pr) - 2.2114
        elif 1.2 <= T_pr <= 1.4:
            dZr = (0.0522 * T_pr) - 0.8511
        elif 1.4 <= T_pr <= 2.0:
            dZr = (0.1391 * T_pr) - 0.2988
        elif 2.0 <= T_pr <= 3.0:
            dZr = (0.0295 * T_pr) - 0.0825
    elif 1.2 <= P_pr <= 2.8:
        if 1.05 <= T_pr <= 1.2:
            dZr = -(1.3570 * T_pr) + 1.4942
        elif 1.2 <= T_pr <= 1.4:
            dZr = (0.1717 * T_pr) - 0.3232
        elif 1.4 <= T_pr <= 2.0:
            dZr = (0.0984 * T_pr) - 0.2053
        elif 2.0 <= T_pr <= 3.0:
            dZr = (0.0211 * T_pr) - 0.0527
    elif 2.8 <= P_pr <= 5.4:
        if 1.05 <= T_pr <= 1.2:
            dZr = -(0.3278 * T_pr) + 0.4752
        elif 1.2 <= T_pr <= 1.4:
            dZr = -(0.2521 * T_pr) + 0.3871
        elif 1.4 <= T_pr <= 2.0:
            dZr = -(0.0284 * T_pr) + 0.0625
        elif 2.0 <= T_pr <= 3.0:
            dZr = -(0.0041 * T_pr) + 0.0039
    elif 5.4 <= P_pr <= 15:
        if 1.05 <= T_pr <= 3.0:
            dZr = (0.711 + (3.66*T_pr))**(-1.4667)
    
    cr = (1/P_pr) - ((1/z)*dZr)
    
    cg = cr/Ppc
    
    return cg




