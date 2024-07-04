"""
### oil_density.py
#### Functions:
    - vazquez_beggs
    - petrosky_farshad
    - ahmed
    - standing
    - rho_o_basic
    
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
def vazquez_beggs(P: float|int, T: float|int, Pb: float|int, Bo: float, Rs: float, gamma_gas: float, gamma_api: float|int) -> float:
    """This is a correlation for determining the oil density, through Vazquez & Beggs.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        Pb (float | int): Bubble Pressure [psia]
        Bo (float): Factor Volumetric Oil [RB/STB]
        Rs (float): Solution Gas-Oil Ratio [PCN/BF]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity

    #### Returns:
        rho_o (float): Oil Density [lb/ft^3]
    """    
    gamma_oil = 141.5 / (gamma_api + 131.5)
        
    rho_ob = ((350*gamma_oil) + (0.0764*gamma_gas*Rs)) / (5.615*Bo)
    
    if P > Pb:
        from math import exp, log
        A = (1e-5)*((-1433)+(5*Rs)+(17.2*(T-460))-(1180*gamma_gas)+(12.61*gamma_api))
        rho_o = rho_ob*exp((-A)*log(P/Pb))
        return round(rho_o, 5)
    
    return round(rho_ob, 5)

@vectorize_decorator
def petrosky_farshad(P: float|int, T: float|int, Pb: float|int, Bo: float, Rs: float, gamma_gas: float, gamma_api: float|int) -> float:
    """This is a correlation for determining the oil density, through Vazquez & Beggs.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        Pb (float | int): Bubble Pressure [psia]
        Bo (float): Factor Volumetric Oil [RB/STB]
        Rs (float): Solution Gas-Oil Ratio [PCN/BF]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity

    #### Returns:
        rho_o (float): Oil Density [lb/ft^3]
    """    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    
    rho_ob = ((350*gamma_oil) + (0.0764*gamma_gas*Rs)) / (5.615*Bo)
    
    if P > Pb:
        from math import exp
        A = (4.1646e-7*(Rs**0.69357))*(gamma_gas**0.1885)*(gamma_api**0.3272)*((T-460)**0.6729)
        rho_o = rho_ob * exp(A*((P**0.4094)-(Pb**0.4094)))
        return round(rho_o, 5)
    
    return round(rho_ob, 5)

@vectorize_decorator
def ahmed(P: float|int, Pb: float|int, Bo: float, Rs: float, gamma_gas: float, gamma_api: float|int) -> float:
    """This is a correlation for determining the oil density, through Ahmed.

    #### Args:
        P (float | int): Pressure [psia]
        Pb (float | int): Bubble Pressure [psia]
        Bo (float): Factor Volumetric Oil [RB/STB]
        Rs (float): Solution Gas-Oil Ratio [PCN/BF]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity

    #### Returns:
        rho_o (float): Oil Density [lb/ft^3]
    """    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    
    rho_ob = ((350*gamma_oil) + (0.0764*gamma_gas*Rs)) / (5.615*Bo)
    
    if P > Pb:
        from math import exp
        B = -(4.588893+(0.025999*Rs))**(-1)
        rho_o = rho_ob*exp(B*(exp(-0.00018473*P)-exp(-0.00018473*Pb)))
        return round(rho_o, 5)
    
    return round(rho_ob, 5)

@vectorize_decorator
def standing(T: float|int, Rs: float, gamma_gas: float, gamma_api: float|int) -> float:
    """This is a correlation for determining the oil density, through Standing.

    #### Args:
        T (float | int): Temperature [oR]
        Rs (float): Solution Gas-Oil Ratio [PCN/BF]
        gamma_gas (float): Specific Gravity Gas
        gamma_api (float | int): API Gravity

    #### Returns:
        rho_o (float): Oil Density [lb/ft^3]
    """    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    
    rho_o = ((62.4*gamma_oil)+(0.0136*Rs*gamma_gas)) / (0.972+(0.000147*(((Rs*((gamma_gas/gamma_oil)**0.5))+(1.25*(T-460)))**1.175)))
    
    return round(rho_o, 5)

@vectorize_decorator
def rho_o_basic(P: float|int, Pb: float|int, Bo: float, Rs: float, gamma_api: float|int, co: float) -> float:
    """The following equation is obtained for integrating of the compressibility definition.

    .. math::
        rho_o = rho_ob*exp(co*(P-Pb))
    
    #### Args:
        P (float | int): Pressure [psia]
        Pb (float | int): Bubble Pressure [psia]
        Bo (float): Factor Volumetric Oil [RB/STB]
        Rs (float): Solution Gas-Oil Ratio [PCN/BF]
        gamma_api (float | int): API Gravity
        co (float): Oil Compressibility [psia^-1]

    #### Returns:
        rho_o (float): Oil Density [lb/ft^3]
    """    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    
    gamma_gd = ((12.5+gamma_api)/50) - (3.5715e-6*gamma_api*Rs)
    
    rho_ob = ((350*gamma_oil) + (0.0764*gamma_gd*Rs)) / (5.615*Bo)
    
    if P > Pb:
        from math import exp
        rho_o = rho_ob * exp(co*(Pb-P))
        return round(rho_o, 5)
    
    return round(rho_ob, 5)


# import matplotlib.pyplot as plt
# from dissolved_gas_oil_ratio import standing as rs_st
# from oil_volumetric_factor import standing as bo_st

# presiones = list(range(500, 3100, 100))
# Pburbuja = 2500
# temperatura = 180 + 460
# api_oil = 31
# grav_gas = 0.95
# compresibilidad = 9.61e-6

# disolubilidad = [rs_st(p, temperatura, grav_gas, api_oil, Pb=Pburbuja) for p in presiones]

# factor_oil = [bo_st(temperatura, rs, api_oil, grav_gas, P=p, Pb=Pburbuja, co=compresibilidad) for rs, p in zip(disolubilidad, presiones)]

# #densidad = [ahmed(p, temperatura, Pburbuja, bo, rs, grav_gas, api_oil) for p, bo, rs in zip(presiones, factor_oil, disolubilidad)]
# #densidad = [rho_oil(p, Pburbuja, bo, rs, api_oil, compresibilidad) for p, bo, rs in zip(presiones, factor_oil, disolubilidad)]
# densidad = [standing(temperatura, rs, grav_gas, api_oil) for rs in disolubilidad]


# # print(presiones)
# # print(disolubilidad)
# # print(factor_oil)
# print(densidad)

# plt.figure(figsize=(10, 6))
# plt.plot(presiones, densidad, marker='o', linestyle='--')

# plt.show()
# plt.grid(True)