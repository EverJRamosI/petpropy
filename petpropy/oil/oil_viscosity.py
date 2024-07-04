"""
### oil_viscosity.py
#### Functions:
    - **Dead Oil Viscosity**
        - beal_muod
        - beggs_robinson_muod
        - glaso_muod
        - egbogad_muod
        - kartoatmodjo_schimdt_muod
    - **Saturated Oil Viscosity**
        - chew_connally_muob
        - beggs_robinson_muob
        - kartoatmodjo_schimdt_muob
    - **Unsaturated Oil Viscosity**
        - beal_muo
        - vazquez_beggs_muo
        - kartoatmodjo_schmidt_muo
    - **mu_oil**: This function is a global of all functions.
    
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
def beal_muod(T: float|int, gamma_api: float|int) -> float:
    """This is a equation for determining the dead oil viscosity, through Beal C.

    #### Args:
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity

    #### Returns:
        mu_od (float): Dead Oil Viscosity [cp]
    """    
    a = 10**(0.43 + (8.33/gamma_api))
    mu_od = (0.32 + ((1.8e7)/(gamma_api**4.53))) * ((360/((T-460)+200))**a)
    return round(mu_od, 5)

@vectorize_decorator
def beggs_robinson_muod(T: float|int, gamma_api: float|int) -> float:
    """This is a equation for determining the dead oil viscosity, through Beggs H.D. & Robinson J.R.

    #### Args:
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity

    #### Returns:
        mu_od (float): Dead Oil Viscosity [cp]
    """    
    z = 3.0324 - (0.02023*gamma_api)
    y = 10**z
    x = y*((T-460)**(-1.163))
    mu_od = (10**x) - 1
    
    return round(mu_od, 5)

@vectorize_decorator
def glaso_muod(T: float|int, gamma_api: float|int) -> float:
    """This is a equation for determining the dead oil viscosity, through Glaso O.

    #### Args:
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity

    #### Returns:
        mu_od (float): Dead Oil Viscosity [cp]
    """    
    from math import log10
    
    mu_od = 3.141e10 * ((T-460)**(-3.444)) * ((log10(gamma_api))**((10.313*log10(T-460))-36.447))
    
    return round(mu_od, 5)

@vectorize_decorator
def egbogad_muod(T: float|int, gamma_api: float|int) -> float:
    """This is a equation for determining the dead oil viscosity, through Egbogah E.O.

    #### Args:
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity

    #### Returns:
        mu_od (float): Dead Oil Viscosity [cp]
    """    
    from math import log10
    
    A = 1.8653 - (0.025086*gamma_api) - (0.5644*log10(T-460))
    B = 10**A
    C = 10**B
    mu_od = C -1
    
    return round(mu_od, 5)

@vectorize_decorator
def kartoatmodjo_schmidt_muod(T: float|int, gamma_api: float|int) -> float:
    """This is a equation for determining the dead oil viscosity, through Kartoatmodjo T. & Schmidt Z.

    #### Args:
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity

    #### Returns:
        mu_od (float): Dead Oil Viscosity [cp]
    """    
    from math import log10
    mu_od = 16.0e8 * ((T-460)**(-2.8177)) * ((log10(gamma_api))**((5.7526*log10(T-460))-26.9718))
    
    return round(mu_od, 5)

@vectorize_decorator
def chew_connally_muob(Rs: float, mu_od: float) -> float:
    """This is a equation for determining the oil viscosity at the bubble point, through Chew J.N. & Connaly C.A.

    #### Args:
        Rs (float): Solution Gas-Oil Ratio [PCN/BF]
        mu_od (float): Dead Oil Viscosity [cp]

    #### Returns:
        mu_ob (float): Oil Viscosity at the Bubble Point [cp]
    """    
    A = 10**(Rs * ((2.2e-7*Rs)-7.4e-4))
    b = (0.68/(10**(8.62e-5*Rs))) + (0.25/(10**(1.1e-3*Rs))) + (0.062/(10**(3.74e-3*Rs)))
    mu_ob = A * (mu_od**b)
    
    return round(mu_ob, 5)

@vectorize_decorator
def beggs_robinson_muob(Rs: float, mu_od: float) -> float:
    """This is a equation for determining the oil viscosity at the bubble point, through Beggs H.D. & Robinson J.R.

    #### Args:
        Rs (float): Solution Gas-Oil Ratio [PCN/BF]
        mu_od (float): Dead Oil Viscosity [cp]

    #### Returns:
        mu_ob (float): Oil Viscosity at the Bubble Point [cp]
    """  
    a = 10.715 * ((Rs + 100)**(-0.515))
    b = 5.44 * ((Rs + 150)**(-0.338))
    mu_ob = a * (mu_od**b)
    
    return round(mu_ob, 5)

@vectorize_decorator
def kartoatmodjo_schmidt_muob(Rs: float, mu_od: float) -> float:
    """This is a equation for determining the oil viscosity at the bubble point, through Kartoatmodjo T. & Schimdt Z.
    
    #### Args:
        Rs (float): Solution Gas-Oil Ratio [PCN/BF]
        mu_od (float): Dead Oil Viscosity [cp]

    #### Returns:
        mu_ob (float): Oil Viscosity at the Bubble Point [cp]
    """  
    b = 10**((-0.00081)*Rs)
    A = (0.2001 + (0.8428 * (10**((-0.000845)*Rs)))) * (mu_od**(0.43 + (0.5165*b)))
    mu_ob = (-0.06821) + (0.9824*A) + (40.34e-5*(A**2))
    
    return round(mu_ob, 5)

@vectorize_decorator
def beal_muo(P: float|int, Pb: float|int, mu_ob: float) -> float:
    """This is a equation for determining the oil viscosity over the bubble point, through Beal C.

    #### Args:
        P (float | int): Pressure [psia]
        Pb (float | int): Bubble Pressure [psia]
        mu_ob (float): Oil Viscosity at the Bubble Point [cp]

    #### Returns:
        mu_o (float): Oil Viscosity [cp]
    """    
    mu_o = (((0.024*(mu_ob**1.6))+(0.038*(mu_ob**0.56)))*(0.001*(P-Pb))) + mu_ob
    
    return round(mu_o, 5)

@vectorize_decorator
def vazquez_beggs_muo(P: float|int, Pb: float|int, mu_ob: float) -> float:
    """This is a equation for determining the oil viscosity over the bubble point, through Vazquez M.E. & Beggs H.D.

    #### Args:
        P (float | int): Pressure [psia]
        Pb (float | int): Bubble Pressure [psia]
        mu_ob (float): Oil Viscosity at the Bubble Point [cp]

    #### Returns:
        mu_o (float): Oil Viscosity [cp]
    """ 
    from math import exp
    
    m = 2.6*(P**1.187)*exp((-11.513)-(8.98e-5*P))
    mu_o = mu_ob * ((P/Pb)**m)
    
    return round(mu_o, 5)

@vectorize_decorator
def kartoatmodjo_schmidt_muo(P: float|int, Pb: float|int, mu_ob: float) -> float:
    """This is a equation for determining the oil viscosity over the bubble point, through Kartoatmodjo T. & Schmidt Z.

    #### Args:
        P (float | int): Pressure [psia]
        Pb (float | int): Bubble Pressure [psia]
        mu_ob (float): Oil Viscosity at the Bubble Point [cp]

    #### Returns:
        mu_o (float): Oil Viscosity [cp]
    """ 
    mu_o = (1.00081*mu_ob) + (1.127e-3*(P-Pb)*((-65.17e-4*(mu_ob**1.8148))+(0.038*(mu_ob**1.59))))
    
    return round(mu_o, 5)

@vectorize_decorator
def mu_oil(P: float|int, T: float|int, gamma_api: float|int, Rs: float, *, Pb: float|int, model_uod: str='Karto', model_uob: str='Karto', model_uo: str='Karto') -> float:
    """This is a function for determining the oil viscosity, it is recommending to use the function because of all functions the oil viscosity it is here.\n
    The following models is:
    
    mu_od:
        - 'Beal': beal_muod
        - 'Beggs': beggs_robinson_muod
        - 'Glaso': glaso_muod
        - 'Egbogad': egbogad_muod
        - 'Karto': kartoatmodjo_schmidt_muod
    
    mu_ob:
        - 'Beggs': beggs_robinson_muob
        - 'Karto': kartoatmodjo_schmidt_muob
        - 'Chew': chew_connally_muob
    
    mu_o:
        - 'Beal': beal_muo
        - 'Karto': kartoatmodjo_schmidt_muo
        - 'Vazquez': vazquez_beggs_muo

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        Rs (float): Solution Gas-Oil Ratio [PCN/BF]
        Pb (float | int): Bubble Pressure [psia]
        model_uod (str, optional): Model of the Dead Oil Viscosity. Defaults to 'Karto'.
        model_uob (str, optional): Model of the Oil Viscosity at the Bubble Point. Defaults to 'Karto'.
        model_uo (str, optional): Oil Viscosity over Bubble Point. Defaults to 'Karto'.

    #### Returns:
        mu_o (float): Oil Viscosity [cp]
    """    
        
    mu_od = {
        'Beal': beal_muod, 
        'Beggs': beggs_robinson_muod, 
        'Glaso': glaso_muod, 
        'Egbogad': egbogad_muod, 
        'Karto': kartoatmodjo_schmidt_muod
        }
    
    mu_ob = {
        'Beggs': beggs_robinson_muob, 
        'Karto': kartoatmodjo_schmidt_muob, 
        'Chew': chew_connally_muob
        }
    mu_o = {
        'Beal': beal_muo, 
        'Karto': kartoatmodjo_schmidt_muo, 
        'Vazquez': vazquez_beggs_muo
        }
    
    if model_uod in mu_od:
        mu_od_value = mu_od[model_uod](T, gamma_api)
        
    if model_uob in mu_ob:
        mu_ob_value = mu_ob[model_uob](Rs, mu_od_value)
        
    if model_uo in mu_o:
        mu_o_value = mu_o[model_uo](P, Pb, mu_ob_value)
        
    return round(mu_o_value, 5)
        

# import matplotlib.pyplot as plt        
# from dissolved_gas_oil_ratio import standing

# presiones = list(range(100, 3100, 100))
# temperatura = 180 + 460
# grav_gas = 0.95
# api_oil = 31

# disolubilidad = [standing(p, temperatura, grav_gas, api_oil, Pb=2500) for p in presiones]

# viscosidades = [mu_oil(p, temperatura, api_oil, r, Pb=2500, model_uod='Beal', model_uo='Beal') for p, r in zip(presiones, disolubilidad)]

# print(presiones)
# print(disolubilidad)
# print(viscosidades)

# plt.figure(figsize=(10, 6))
# plt.plot(presiones, viscosidades, marker='o', linestyle='--')

# plt.show()
# plt.grid(True)

#print('beal', beal_muod((180+460), 31))
#print('beggs', beggs_robinson_muod((180+460), 31))
#print('glaso', glaso_muod((180+460), 31))
#print('egboad', egbogad_muod((180+460), 31))
#print('karto', kartoatmodjo_schmidt_muod((180+460), 31))

#print('beggs', beggs_robinson_muob(675, 2.65))
#print('karto', kartoatmodjo_schmidt_muob(675, 2.65))
#print('chew', chew_connally_muob(675, 2.65))

#print('beal', beal_muo(4000, 2500, 0.74))
#print('karto', kartoatmodjo_schmidt_muo(4000, 2500, 0.74))
#print('vazquez', vazquez_beggs_muo(4000, 2500, 0.74))

