"""
### solution_gas_oil_ratio.py
#### Functions:
    - standing
    - lasater
    - vazquez
    - glaso
    - total
    - almarhoun
    - dokla_osman
    - petrosky_farshad
    - kartoatmodjo_schmidt
    
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
def standing(P: float|int, T: float|int, gamma_gas: float, gamma_api: float|int, *, Pb: float|int) -> float:
    """This is a correlation for determining Solution Gas-Oil Ratio, through Standing M.B.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity
        Pb (float | int): Bubble Pressure [oR]

    #### Returns:
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]
    """    
    if P >= Pb:
        P = Pb
        Rs = gamma_gas * ((((P/18.2)+1.4)*(10**((0.0125*gamma_api)-(0.00091*(T-460)))))**1.2048)
        
        return Rs
    
    if P < Pb:
        Rs = gamma_gas * ((((P/18.2)+1.4)*(10**((0.0125*gamma_api)-(0.00091*(T-460)))))**1.2048)
    
        return Rs

@vectorize_decorator
def lasater(P: float|int, T: float|int, gamma_gas: float, gamma_api: float|int, *, Pb: float|int) -> float:
    """This is a correlation for determining Solution Gas-Oil Ratio, through Lasater J.A.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity
        Pb (float | int): Bubble Pressure [oR]

    #### Returns:
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]
    """  
    from math import log
    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    
    if P >= Pb:
        P = Pb
        Pf = (Pb*gamma_gas)/T
        if Pf < 3.29:
            y_gas = 0.359 * log(((1.473*P*gamma_gas)/T)+0.476)
        elif Pf >= 3.29:
            y_gas = (((0.121*P*gamma_gas)/T)-0.236)**0.281
        
        if gamma_api <= 40:
            m_oil = 630 - (10 * gamma_api)
        elif gamma_api > 40:
            m_oil = 73110 / (gamma_api**1.562)
        
        Rs = (132755*gamma_oil*y_gas) / (m_oil*(1-y_gas))
    
        return Rs
    
    if P < Pb:
        Pb = P
        Pf = (Pb*gamma_gas)/T
        if Pf < 3.29:
            y_gas = 0.359 * log(((1.473*P*gamma_gas)/T)+0.476)
        elif Pf >= 3.29:
            y_gas = (((0.121*P*gamma_gas)/T)-0.236)**0.281
        
        if gamma_api <= 40:
            m_oil = 630 - (10 * gamma_api)
        elif gamma_api > 40:
            m_oil = 73110 / (gamma_api**1.562)
        Rs = (132755*gamma_oil*y_gas) / (m_oil*(1-y_gas))
    
        return Rs

@vectorize_decorator
def vazquez(P: float|int, T: float|int, gamma_gas: float, gamma_api: float|int, *, Pb: float|int, P_sp: float=0, T_sp: float=0) -> float:
    """This is a correlation for determining Solution Gas-Oil Ratio, through Vazquez M.E. & Beggs H.D.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity
        Pb (float | int): Bubble Pressure [oR]
        P_sp (float, optional): Separator Pressure [psia]. Defaults to 0.
        T_sp (float, optional): Separator Temperature [oR]. Defaults to 0.

    #### Returns:
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]
    """   
    from math import log10, exp
    
    if gamma_api <= 30:
        C1 = 0.0362
        C2 = 1.0937
        C3 = 25.724
    elif gamma_api > 30:
        C1 = 0.0178
        C2 = 1.1870
        C3 = 23.931
        
    if P >= Pb:
        P = Pb
        if P_sp != 0 and T_sp != 0:
            gamma_gc = gamma_gas * (1 + (5.912e-5 * gamma_api * T_sp * log10(P_sp/114.7)))
            Rs = C1 * gamma_gc * (P**C2) * exp((C3*gamma_api)/(T))
        else:
            Rs = C1 * gamma_gas * (P**C2) * exp((C3*gamma_api)/(T))
    
        return Rs
    
    if P < Pb:
        if P_sp != 0 and T_sp != 0:
            gamma_gc = gamma_gas * (1 + (5.912e-5 * gamma_api * T_sp * log10(P_sp/114.7)))
            Rs = C1 * gamma_gc * (P**C2) * exp((C3*gamma_api)/(T))
        else:
            Rs = C1 * gamma_gas * (P**C2) * exp((C3*gamma_api)/(T))
    
        return Rs

@vectorize_decorator
def glaso(P: float|int, T: float|int, gamma_gas: float, gamma_api: float|int, *, Pb: float|int) -> float:
    """This is a correlation for determining Solution Gas-Oil Ratio, through Glaso O.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity
        Pb (float | int): Bubble Pressure [oR]

    #### Returns:
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]
    """
    from math import log10
        
    if P >= Pb:
        P = Pb
        F = 10**(2.8869 - ((14.1811-(3.3093*log10(P)))**0.5))
        Rs = gamma_gas * ((F*((gamma_api**0.989) / ((T-460)**0.172)))**1.2255)
        
        return Rs
    
    if P < Pb:
        F = 10**(2.8869 - ((14.1811-(3.3093*log10(P)))**0.5))
        Rs = gamma_gas * ((F*((gamma_api**0.989) / ((T-460)**0.172)))**1.2255)
    
        return Rs

@vectorize_decorator
def total(P: float|int, T: float|int, gamma_gas: float, gamma_api: float|int, *, Pb: float|int) -> float:
    """This is a correlation for determining Solution Gas-Oil Ratio, through TOTAL C.F.P.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity
        Pb (float | int): Bubble Pressure [oR]

    #### Returns:
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]
    """
    if gamma_api <= 10:
        C1 = 12.2651
        C2 = 0.030405
        C3 = 0
        C4 = 0.9669
    elif 10 < gamma_api <= 35:
        C1 = 15.0057
        C2 = 0.0152
        C3 = 4.484e-4
        C4 = 1.0950
    elif 35 < gamma_api <= 45:
        C1 = 112.925
        C2 = 0.0248
        C3 = -1.469e-3
        C4 = 1.1290
        
    if P >= Pb:
        P = Pb
        Rs = gamma_gas * (((P/C1)*(10**((C2*gamma_api)-(C3*(T-460)))))**C4)
    
        return Rs
    
    if P < Pb:
        Rs = gamma_gas * (((P/C1)*(10**((C2*gamma_api)-(C3*(T-460)))))**C4)
    
        return Rs

@vectorize_decorator
def almarhoun(P: float|int, T: float|int, gamma_gas: float, gamma_api: float|int, *, Pb: float|int) -> float:
    """This is a correlation for determining Solution Gas-Oil Ratio, through Al-Marhoun M.A.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity
        Pb (float | int): Bubble Pressure [oR]

    #### Returns:
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]
    """
    gamma_oil = 141.5 / (gamma_api + 131.5)
    
    if P >= Pb:
        P = Pb
        Rs = (185.84321 * P * (gamma_gas**1.87784) * (gamma_oil**(-3.1437) * (T**(-1.32657))))**1.3984
    
        return Rs
    
    if P < Pb:
        Rs = (185.84321 * P * (gamma_gas**1.87784) * (gamma_oil**(-3.1437) * (T**(-1.32657))))**1.3984
    
        return Rs

@vectorize_decorator
def dokla_osman(P: float|int, T: float|int, gamma_gas: float, gamma_api: float|int, *, Pb: float|int) -> float:
    """This is a correlation for determining Solution Gas-Oil Ratio, through Dokla M.E. & Osman M.E.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity
        Pb (float | int): Bubble Pressure [oR]

    #### Returns:
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]
    """
    gamma_oil = 141.5 / (gamma_api + 131.5)
    
    if P >= Pb:
        P = Pb
        Rs = (0.11956e-3 * P * (gamma_gas**1.01049) * (gamma_oil**(-0.107991)) * (T**0.952584))**1.3811
    
        return Rs
    
    if P < Pb:
        Rs = (0.11956e-3 * P * (gamma_gas**1.01049) * (gamma_oil**(-0.107991)) * (T**0.952584))**1.3811
    
        return Rs

@vectorize_decorator
def petrosky_farshad(P: float|int, T: float|int, gamma_gas: float, gamma_api: float|int, *, Pb: float|int) -> float:
    """This is a correlation for determining Solution Gas-Oil Ratio, through Petrosky G.E. & Farshad F.F.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity
        Pb (float | int): Bubble Pressure [oR]

    #### Returns:
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]
    """
    if P >= Pb:
        P = Pb
        Rs = ((gamma_gas**0.8439) * ((P/112.727)+12.34) * (10**((7.916e-4*(gamma_api**1.5410))-(4.561e-5*((T-460)**1.3911)))))**1.73184
    
        return Rs
    
    if P < Pb:
        Rs = ((gamma_gas**0.8439) * ((P/112.727)+12.34) * (10**((7.916e-4*(gamma_api**1.5410))-(4.561e-5*((T-460)**1.3911)))))**1.73184
    
        return Rs

@vectorize_decorator
def kartoatmodjo_schmidt(P: float|int, T: float|int, gamma_gas: float, gamma_api: float|int, *, Pb: float|int, P_sp: float=0, T_sp: float=0) -> float:
    """This is a correlation for determining Solution Gas-Oil Ratio, through Kartoatmodjo T. & Schmidt Z.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity
        Pb (float | int): Bubble Pressure [oR]
        P_sp (float, optional): Separator Pressure [psia]. Defaults to 0.
        T_sp (float, optional): Separator Temperature [oR]. Defaults to 0.

    #### Returns:
        Rs (float): Solution Gas-Oil Ratio [PCN/BN]
    """
    from math import log10
            
    if gamma_api <= 30:
        C1 = 0.05958
        C2 = 0.7972
        C3 = 13.1405
        C4 = 0.9986
    elif gamma_api > 30:
        C1 = 0.03150
        C2 = 0.7587
        C3 = 11.2895
        C4 = 0.9143
        
    if P >= Pb:
        P = Pb
        if P_sp != 0 and T_sp != 0:
            gamma_gc = gamma_gas * (1 + (0.1595 * (gamma_api**0.4078) * (T_sp**(-0.2466)) * log10(P_sp/114.7)))
            Rs = C1 * (gamma_gc**C2) * (P**(1/C4)) * (10**((C3*gamma_api)/(T)))
        else:
            Rs = C1 * (gamma_gas**C2) * (P**(1/C4)) * (10**((C3*gamma_api)/(T)))
    
        return Rs
    
    if P < Pb:
        if P_sp != 0 and T_sp != 0:
            gamma_gc = gamma_gas * (1 + (0.1595 * (gamma_api**0.4078) * (T_sp**(-0.2466)) * log10(P_sp/114.7)))
            Rs = C1 * (gamma_gc**C2) * (P**(1/C4)) * (10**((C3*gamma_api)/(T)))
        else:
            Rs = C1 * (gamma_gas**C2) * (P**(1/C4)) * (10**((C3*gamma_api)/(T)))
    
        return Rs




#print(rs_oil(3000, 460, 0.65, 31, Pb=2500, model='Standing'))
# print("Standing", standing(4000, (180+460), 0.95, 31, Pb=2500))
# print('Lasater', lasater(4000, (180+460), 0.95, 31, Pb=2500))
# print('Vazquez', vazquez(4000, (180+460), 0.95, 31, Pb=2500))
# print('Glaso', glaso(4000, (180+460), 0.95, 31, Pb=2500))
# print('total', total(4000, (180+460), 0.95, 31, Pb=2500))
# print('almarhoun', almarhoun(4000, (180+460), 0.95, 31, Pb=2500))
# print('Dokla', dokla_osman(4000, (180+460), 0.95, 31, Pb=2500))
# print('petrosky', petrosky_farshad(4000, (180+460), 0.95, 31, Pb=2500))
# print('kartoatmodjo', kartoatmodjo_schmidt(4000, (180+460), 0.95, 31, Pb=2500))
# print("--------------")
# print("Standing", standing(2000, (180+460), 0.95, 31, Pb=2500))
# print('Lasater', lasater(2000, (180+460), 0.95, 31, Pb=2500))
# print('Vazquez', vazquez(2000, (180+460), 0.95, 31, Pb=2500))
# print('Glaso', glaso(2000, (180+460), 0.95, 31, Pb=2500))
# print('total', total(2000, (180+460), 0.95, 31, Pb=2500))
# print('almarhoun', almarhoun(2000, (180+460), 0.95, 31, Pb=2500))
# print('Dokla', dokla_osman(2000, (180+460), 0.95, 31, Pb=2500))
# print('petrosky', petrosky_farshad(2000, (180+460), 0.95, 31, Pb=2500))
# print('kartoatmodjo', kartoatmodjo_schmidt(2000, (180+460), 0.95, 31, Pb=2500))

#import numpy as np
#import matplotlib.pyplot as plt
#
#presiones = np.arange(500, 3100, 100)
#temperatura = 180 + 460
#gam_esp = 0.95
#api_oil = 31
#pb = 2500
#compr = 9.61e-6
#
#rs = [standing(p, temperatura, gam_esp, api_oil, Pb=pb) for p in presiones]
#rs = [lasater(p, temperatura, gam_esp, api_oil, Pb=pb) for p in presiones]
#rs = [vazquez(p, temperatura, gam_esp, api_oil, Pb=pb) for p in presiones]
#rs = [glaso(p, temperatura, gam_esp, api_oil, Pb=pb) for p in presiones]
#rs = [total(p, temperatura, gam_esp, api_oil, Pb=pb) for p in presiones]
#rs = [almarhoun(p, temperatura, gam_esp, api_oil, Pb=pb) for p in presiones]
#rs = [dokla_osman(p, temperatura, gam_esp, api_oil, Pb=pb) for p in presiones]
#rs = [petrosky_farshad(p, temperatura, gam_esp, api_oil, Pb=pb) for p in presiones]
#rs = [kartoatmodjo_schmidt(p, temperatura, gam_esp, api_oil, Pb=pb) for p in presiones]
#
#
#print(presiones)
#print(rs)
#
#
#plt.figure(figsize=(10, 6))
#plt.plot(presiones, rs, marker='o', linestyle='--')
#
#plt.show()
#plt.grid(True)
