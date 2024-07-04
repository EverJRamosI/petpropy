"""
### oil_volumetric_factor.py
#### Functions:
    - standing
    - vazquez_beggs
    - glaso
    - total
    - almarhoun
    - dokla_osman
    - petrosky_farshad
    - kartoatmodjo_schmidt
    - factor_oil
    
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
def standing(T: float|int, R_sb: float, gamma_api: float|int, gamma_gas: float, *, P: float=0, Pb: float=0, co: float=0) -> float:
    """This is a correlation for determining the oil volumetric factor, through Standing M. B.

    #### Args:
        T (float | int): Temperature [oR]
        R_sb (float): Solution Gas-Oil Ratio [PCN/BN]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        P (float, optional): Pressure [psia]. Defaults to 0.
        Pb (float, optional): Bubble Pressure [psia]. Defaults to 0.
        co (float, optional): Oil Compressibility [psia^-1]. Defaults to 0.

    #### Returns:
        Bo (float): Oil Volumetric Factor [RB/STB]
    """    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    F = (R_sb * ((gamma_gas/gamma_oil)**(0.5))) + (1.25 * (T-460))
    B_ob = 0.9759 + (12e-5 * (F**1.2))
    
    if P != 0 and Pb != 0 and co != 0:
        from math import exp
        Bo = B_ob*exp(co*(Pb-P))
        
        return Bo
    
    return B_ob

@vectorize_decorator
def vazquez_beggs(T: float|int, R_sb: float, gamma_api: float|int, gamma_gas: float, *, P: float=0, Pb: float=0, co: float=0, P_sp: float=0, T_sp: float=0) -> float:
    """This is a correlation for determining the oil volumetric factor, through Vazquez M.E. & Beggs H.D.

    #### Args:
        T (float | int): Temperature [oR]
        R_sb (float): Solution Gas-Oil Ratio [PCN/BN]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        P (float, optional): Pressure [psia]. Defaults to 0.
        Pb (float, optional): Bubble Pressure [psia]. Defaults to 0.
        co (float, optional): Oil Compressibility [psia^-1]. Defaults to 0.
        P_sp (float, optional): Separator Pressure [psia]. Defaults to 0.
        T_sp (float, optional): Separator Temperature [oR]. Defaults to 0.

    #### Returns:
        Bo (float): Oil Volumetric Factor [RB/STB]
    """    
    from math import log10, exp
        
    if gamma_api <= 30:
        C1 = 4.677e-4
        C2 = 1.751e-5
        C3 = -1.8106e-8
    elif gamma_api > 30:
        C1 = 4.670e-4
        C2 = 1.100e-5
        C3 = 1.3370e-9
    
    if P_sp != 0 and T_sp != 0:
        T_sp -= 460
        gamma_gc = gamma_gas * (1 + (0.1595 * (gamma_api**0.4078) * (T_sp**(-0.2466)) * log10(P_sp/114.7)))   
        B_ob = 1 + (C1*R_sb) + (C2*((T-460)-60)*(gamma_api/gamma_gc)) + (C3*R_sb*((T-460)-60)*(gamma_api/gamma_gc))
    else:
        B_ob = 1 + (C1*R_sb) + (C2*((T-460)-60)*(gamma_api/gamma_gas)) + (C3*R_sb*((T-460)-60)*(gamma_api/gamma_gas))
    
    if P != 0 and Pb != 0 and co != 0:
        Bo = B_ob*exp(co*(Pb-P))
        
        return Bo
    
    return B_ob

@vectorize_decorator
def glaso(T: float|int, R_sb: float, gamma_api: float|int, gamma_gas: float, *, P: float=0, Pb: float=0, co: float=0) -> float:
    """This is a correlation for determining the oil volumetric factor, through Glaso O.

    #### Args:
        T (float | int): Temperature [oR]
        R_sb (float): Solution Gas-Oil Ratio [PCN/BN]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        P (float, optional): Pressure [psia]. Defaults to 0.
        Pb (float, optional): Bubble Pressure [psia]. Defaults to 0.
        co (float, optional): Oil Compressibility [psia^-1]. Defaults to 0.

    #### Returns:
        Bo (float): Oil Volumetric Factor [RB/STB]
    """    
    from math import log10, exp
    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    F = (R_sb*((gamma_gas/gamma_oil)**0.526)) + (0.968*(T-460))
    B_ob = 1 + (10**((-6.58511)+(2.91329*log10(F))-(0.27683*((log10(F))**2))))
    
    if P != 0 and Pb != 0 and co != 0:
        Bo = B_ob*exp(co*(Pb-P))
        
        return Bo
    
    return B_ob

@vectorize_decorator
def total(T: float|int, R_sb: float, gamma_api: float|int, gamma_gas: float, *, P: float=0, Pb: float=0, co: float=0) -> float:
    """This is a correlation for determining the oil volumetric factor, through TOTAL C.F.P.

    #### Args:
        T (float | int): Temperature [oR]
        R_sb (float): Solution Gas-Oil Ratio [PCN/BN]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        P (float, optional): Pressure [psia]. Defaults to 0.
        Pb (float, optional): Bubble Pressure [psia]. Defaults to 0.
        co (float, optional): Oil Compressibility [psia^-1]. Defaults to 0.

    #### Returns:
        Bo (float): Oil Volumetric Factor [RB/STB]
    """    
    B_ob = 1.022 + (4.857e-4*R_sb) - (2.009e-6*((T-460)-60)*(gamma_api/gamma_gas)) + (17.569e-9*R_sb*((T-460)-60)*(gamma_api/gamma_gas))
    
    if P != 0 and Pb != 0 and co != 0:
        from math import exp
        Bo = B_ob*exp(co*(Pb-P))
        
        return Bo
    
    return B_ob

@vectorize_decorator
def almarhoun(T: float|int, R_sb: float, gamma_api: float|int, gamma_gas: float, *, P: float=0, Pb: float=0, co: float=0) -> float:
    """This is a correlation for determining the oil volumetric factor, through Al-Marhoun M.A.

    #### Args:
        T (float | int): Temperature [oR]
        R_sb (float): Solution Gas-Oil Ratio [PCN/BN]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        P (float, optional): Pressure [psia]. Defaults to 0.
        Pb (float, optional): Bubble Pressure [psia]. Defaults to 0.
        co (float, optional): Oil Compressibility [psia^-1]. Defaults to 0.

    #### Returns:
        Bo (float): Oil Volumetric Factor [RB/STB]
    """    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    F = (R_sb**0.74239) * (gamma_gas**0.323294) * (gamma_oil**(-1.20204))
    B_ob = 0.497069 + (0.862963e-3*T) + (0.182594e-2*F) + (0.318099e-5*(F**2))
    
    if P != 0 and Pb != 0 and co != 0:
        from math import exp
        Bo = B_ob*exp(co*(Pb-P))
        
        return Bo
    
    return B_ob

@vectorize_decorator
def dokla_osman(T: float|int, R_sb: float, gamma_api: float|int, gamma_gas: float, *, P: float=0, Pb: float=0, co: float=0) -> float:
    """This is a correlation for determining the oil volumetric factor, through Dokla M.E. & Osman M.E.

    #### Args:
        T (float | int): Temperature [oR]
        R_sb (float): Solution Gas-Oil Ratio [PCN/BN]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        P (float, optional): Pressure [psia]. Defaults to 0.
        Pb (float, optional): Bubble Pressure [psia]. Defaults to 0.
        co (float, optional): Oil Compressibility [psia^-1]. Defaults to 0.

    #### Returns:
        Bo (float): Oil Volumetric Factor [RB/STB]
    """    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    F = (R_sb**0.773572) * (gamma_gas**0.40402) * (gamma_oil**(-0.882605))
    B_ob = (0.431936e-1) + (0.156667e-2*T) + (0.139775e-2*F) + (0.380525e-5*(F**2))
    
    if P != 0 and Pb != 0 and co != 0:
        from math import exp
        Bo = B_ob*exp(co*(Pb-P))
        
        return Bo
    
    return B_ob

@vectorize_decorator
def petrosky_farshad(T: float|int, R_sb: float, gamma_api: float|int, gamma_gas: float, *, P: float=0, Pb: float=0, co: float=0) -> float:
    """This is a correlation for determining the oil volumetric factor, through Petrosky G.E. & Farshad F.F.

    #### Args:
        T (float | int): Temperature [oR]
        R_sb (float): Solution Gas-Oil Ratio [PCN/BN]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        P (float, optional): Pressure [psia]. Defaults to 0.
        Pb (float, optional): Bubble Pressure [psia]. Defaults to 0.
        co (float, optional): Oil Compressibility [psia^-1]. Defaults to 0.

    #### Returns:
        Bo (float): Oil Volumetric Factor [RB/STB]
    """    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    F = (R_sb**0.338) * ((gamma_gas**0.2914)/(gamma_oil**0.6265)) + (0.24626*(T**0.5371))
    B_ob = 1.0113 + (7.2046e-5*(F**3.0936))
    
    if P != 0 and Pb != 0 and co != 0:
        from math import exp
        Bo = B_ob*exp(co*(Pb-P))
        
        return Bo
    
    return B_ob

@vectorize_decorator
def kartoatmodjo_schmidt(T: float|int, R_sb: float, gamma_api: float|int, gamma_gas: float, *, P: float=0, Pb: float=0, co: float=0, P_sp: float=0, T_sp: float=0) -> float:
    """This is a correlation for determining the oil volumetric factor, through Kartoatmodjo T. & Schmidt Z.

    #### Args:
        T (float | int): Temperature [oR]
        R_sb (float): Solution Gas-Oil Ratio [PCN/BN]
        gamma_api (float | int): API Gravity
        gamma_gas (float): Specific Gravity Gas
        P (float, optional): Pressure [psia]. Defaults to 0.
        Pb (float, optional): Bubble Pressure [psia]. Defaults to 0.
        co (float, optional): Oil Compressibility [psia^-1]. Defaults to 0.
        P_sp (float, optional): Separator Pressure [psia]. Defaults to 0.
        T_sp (float, optional): Separator Temperature [oR]. Defaults to 0.

    #### Returns:
        Bo (float): Oil Volumetric Factor [RB/STB]
    """    
    from math import log10, exp
    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    
    if P_sp != 0 and T_sp != 0:
        gamma_gc = gamma_gas * (1 + (0.1595 * (gamma_api**0.4078) * (T_sp**(-0.2466)) * log10(P_sp/114.7)))
        F = ((R_sb**0.755) * (gamma_gc**0.25) * (gamma_oil**(-1.5))) + (0.45*(T-460))
        B_ob = 0.98496 + (1e-4*(F**1.5))
    else:
        F = ((R_sb**0.755) * (gamma_gas**0.25) * (gamma_oil**(-1.5))) + (0.45*(T-460))
        B_ob = 0.98496 + (1e-4*(F**1.5))
        
    if P != 0 and Pb != 0 and co != 0:
        
        Bo = B_ob*exp(co*(Pb-P))
        
        return Bo
    
    return B_ob

@vectorize_decorator
def factor_oil(P: float|int, Pb: float|int, B_ob: float, co: float) -> float:
    """Equation for determining the oil volumetric factor.

    .. math::
        Bo = Bob*exp(co*(Pb-P))
    #### Args:
        P (float | int): Pressure [psia]
        Pb (float | int): Bubble Pressure [psia]
        B_ob (float): Oil Volumetric Factor at the Bubble Point [RB/STB]
        co (float): Oil Compressibility [psia^-1]

    #### Returns:
        Bo (float): Oil Volumetric Factor [RB/STB]
    """    
    from math import exp
    
    Bo = B_ob*exp(co*(Pb-P))
    
    return Bo




#print('Standing 2500', standing((180+460), 673, 31, 0.95, P=3000, Pb=2500, co=9.61e-6))
#print('Standing 3000', factor_oil(3000, 2500, standing((180+460), 673, 31, 0.95), 9.61e-6))

#import numpy as np
#from dissolved_gas_oil_ratio import standing as st_Rs
#import matplotlib.pyplot as plt
##
#presiones = np.arange(500, 3100, 100)
#temperatura = 180 + 460
#gam_esp = 0.95
#api_oil = 31
#pb = 2500
#compr = 9.61e-6
#
#rs = [st_Rs(p, temperatura, gam_esp, api_oil, Pb=pb) for p in presiones]

#bo = [standing(temperatura, r, api_oil, gam_esp, P=p, Pb=pb, co=compr) for r, p in zip(rs, presiones)]
#bo = [vazquez_beggs(temperatura, r, api_oil, gam_esp, P=p, Pb=pb, co=compr) for r, p in zip(rs, presiones)]
#bo = [glaso(temperatura, r, api_oil, gam_esp, P=p, Pb=pb, co=compr) for r, p in zip(rs, presiones)]
#bo = [total(temperatura, r, api_oil, gam_esp, P=p, Pb=pb, co=compr) for r, p in zip(rs, presiones)]
#bo = [almarhoun(temperatura, r, api_oil, gam_esp, P=p, Pb=pb, co=compr) for r, p in zip(rs, presiones)]
#bo = [dokla_osman(temperatura, r, api_oil, gam_esp, P=p, Pb=pb, co=compr) for r, p in zip(rs, presiones)]
#bo = [petrosky_farshad(temperatura, r, api_oil, gam_esp, P=p, Pb=pb, co=compr) for r, p in zip(rs, presiones)]
#bo = [kartoatmodjo_schmidt(temperatura, r, api_oil, gam_esp, P=p, Pb=pb, co=compr) for r, p in zip(rs, presiones)]
#bo = [b_oil(temperatura, r, api_oil, gam_esp, P=p, Pb=pb, co=compr, model='Total') for r, p in zip(rs, presiones)]

#print(presiones)
#print(rs)
#print(bo)
#
#plt.figure(figsize=(10, 6))
#plt.plot(presiones, bo, marker='o', linestyle='--')
##
#plt.show()
#plt.grid(True)




