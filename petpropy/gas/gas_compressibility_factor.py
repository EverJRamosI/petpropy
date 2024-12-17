"""
### gas_compressibility_factor.py

#### Functions:
    - wichert_aziz_correction
    - wichert_aziz
    - papay
    - brill_beggs
    - hall_yarborough
    - gopal
    - dranchuk_abou_kassem
    - dranchuk_purvis_robinson

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

def wichert_aziz_correction(Ppc, Tpc):
    """This is a correction function of Whichert and Aziz. Corrects the pseudocritical Pressure and Temperature affected by CO2 & H2S.
    
    #### Args:
        - Ppc (float|int): Pressure pseudocritical uncorrected [psia]
        - Tpc (float|int): Temperature pseudocritical uncorrected [oR]

    #### Returns:
        Dict : return Ppc and Tpc corrected in a dictionary
    """    
    y_impurities = ['CO2', 'H2S']
        
    f_impurity = [float(imp.strip()) for imp in input(f'Enter the following impurities {y_impurities} in decimals in the corresponding order and separated by a comma (if you do not have an impurity, enter zero [0]): \nIngrese las siguientes impurezas {y_impurities} en decimales en el orden correspondiente y separados por una coma (si no tiene una impureza coloque cero [0]): ').split(',')]
    
    y_impurities = dict(zip(y_impurities, f_impurity))
    
    varepsilon = 120 * (((y_impurities.get('CO2') + y_impurities.get('H2S'))**0.9) - ((y_impurities.get('CO2') + y_impurities.get('H2S'))**1.6)) + 15 * (((y_impurities.get('H2S'))**0.5) - ((y_impurities.get('H2S'))**4))
    
    Tpc_corr = Tpc - varepsilon
    Ppc_corr = (Ppc*Tpc_corr) / (Tpc + (y_impurities.get('H2S') * (1 - y_impurities.get('H2S')) * varepsilon))
    
    return {"Ppc_corregida=": Ppc_corr, "Tpc_corregida=": Tpc_corr}

def wichert_aziz(Ppc, Tpc, *, yCO2, yH2S):
    """This is a correction function of Whichert and Aziz. Corrects the pseudocritical Pressure and Temperature affected by CO2 & H2S.
    
    #### Args:
        - Ppc (float|int): Pressure pseudocritical uncorrected [psia]
        - Tpc (float|int): Temperature pseudocritical uncorrected [oR]

    #### Returns:
        List : return Ppc and Tpc corrected in a list [Ppc_corr, Tpc_corr]
    """    
    varepsilon = 120*((yCO2+yH2S)**0.9 - (yCO2+yH2S)**1.6) + 15*(yH2S**0.5 - yH2S**4)
    
    Tpc_corr = Tpc - varepsilon
    
    Ppc_corr = (Ppc*Tpc_corr) / (Tpc + (yH2S* (1 - yH2S)) * varepsilon)
    
    return [Ppc_corr, Tpc_corr]

@vectorize_decorator
def papay(P: float, T: float, Ppc: float, Tpc: float) -> float:
    """This is a correlation for determining the z factor for method of Papay. 

    #### Args:
        - P (float): Pressure [psia]
        - T (float): Temperature [oR]
        - Ppc (float): Pressure pseudocritical [psia]
        - Tpc (float): Temperature pseudocritical [oR]

    #### Returns:
        z (float): Gas Compressibility Factor [dimensionless]
    """    
    P_pr = P / Ppc
    T_pr = T / Tpc
    return 1 - ((3.52*P_pr)/(10**(0.9813*T_pr))) + ((0.274*(P_pr**2))/(10**(0.8157*T_pr)))

@vectorize_decorator
def brill_beggs(P: float, T: float, Ppc: float, Tpc: float) -> float:
    """This is a correlation for determining the z factor for method of Brill J.P. & Beggs H.D. 

    #### Args:
        - P (float): Pressure [psia]
        - T (float): Temperature [oR]
        - Ppc (float): Pressure pseudocritical [psia]
        - Tpc (float): Temperature pseudocritical [oR]

    #### Returns:
        z (float): Gas Compressibility Factor [dimensionless]
    """  
    import numpy as np
    
    P_pr = P / Ppc
    T_pr = T / Tpc
    
    A = 1.39 * (T_pr-0.92)**0.5 - 0.36 * T_pr - 0.10
    B = ((0.62-(0.23*T_pr))*P_pr) + (((0.066/(T_pr-0.86))-0.037)*(P_pr**2)) + ((0.32/(10**(9*(T_pr-1))))*(P_pr**6))
    C = 0.132 - (0.32 * np.log10(T_pr))
    D = 10**(0.3106-(0.49*T_pr)+0.1824*(T_pr**2))
    
    with np.errstate(over='ignore'):
        exp_B = np.exp(B)
        
    if np.isinf(exp_B):
        z = A + (C*(P_pr**D))
    else: 
        z = A + ((1-A)/(np.exp(B))) + (C*(P_pr**D))
        
    return z

@vectorize_decorator
def hall_yarborough(P: float, T: float, Ppc: float, Tpc: float) -> float:
    """This is a correlation for determining the z factor for method of Hall K.R. & Yarborough L.

    #### Args:
        - P (float): Pressure [psia]
        - T (float): Temperature [oR]
        - Ppc (float): Pressure pseudocritical [psia]
        - Tpc (float): Temperature pseudocritical [oR]

    #### Returns:
        z (float): Gas Compressibility Factor [dimensionless]
    """  
    import numpy as np
    
    P_pr = P / Ppc
    T_pr = T / Tpc
    t = 1 / T_pr
    
    A = 0.06123 * t * np.exp(-1.2 * ((1 - t)**2))
    B = 14.76 * t - 9.76 * (t**2) + 4.85 * (t**3)
    C = 90.7 * t - 242.2 * (t**2) + 42.4 * (t**3)
    D = 2.18 + 2.82 * t
        
    y1 = 0.00001
    tolerance = 1e-5
    max_iterations = 1000
    iteration = 0
    
    while iteration < max_iterations:
        
        F1 = -A * P_pr + ((y1 + (y1**2) + (y1**3) - y1**4) / ((1 - y1)**3)) - B * (y1**2) + C * (y1**D)
        dF = ((1 + 4 * y1 + 4 * (y1**2) - 4 * (y1**3) + (y1**4)) / ((1 - y1)**4)) - 2 * B * y1 + C * D * (y1**(D - 1))
        
        y1_new = y1 - F1 / dF
        
        y1 = y1_new
        
        if abs(F1) < tolerance:
            break
        
        iteration += 1
        
    z = (0.06125 * P_pr * t * np.exp(-1.2 * ((1 - t) ** 2))) / y1
    
    return z

@vectorize_decorator
def gopal(P: float, T: float, Ppc: float, Tpc: float) -> float:
    """This is a correlation for determining the z factor for method of Gopal V.N.

    #### Args:
        - P (float): Pressure [psia]
        - T (float): Temperature [oR]
        - Ppc (float): Pressure pseudocritical [psia]
        - Tpc (float): Temperature pseudocritical [oR]

    #### Returns:
        z (float): Gas Compressibility Factor [dimensionless]
    """  
    P_pr = P / Ppc
    T_pr = T / Tpc
    
    if 0.2 <= P_pr <= 1.2:
        if 1.05 <= T_pr <= 1.2:
            z = (P_pr * ((1.6643 * T_pr) - 2.2114)) - (0.3647 * T_pr) + 1.4385
        elif 1.2 <= T_pr <= 1.4:
            z = (P_pr * ((0.0522 * T_pr) - 0.8511)) - (0.0364 * T_pr) + 1.0490
        elif 1.4 <= T_pr <= 2.0:
            z = (P_pr * ((0.1391 * T_pr) - 0.2988)) + (0.0007 * T_pr) + 0.9969
        elif 2.0 <= T_pr <= 3.0:
            z = (P_pr * ((0.0295 * T_pr) - 0.0825)) + (0.0009 * T_pr) + 0.9967
    elif 1.2 <= P_pr <= 2.8:
        if 1.05 <= T_pr <= 1.2:
            z = (P_pr * ((-1.3570 * T_pr) + 1.4942)) + (4.6315 * T_pr) - 4.7009
        elif 1.2 <= T_pr <= 1.4:
            z = (P_pr * ((0.1717 * T_pr) - 0.3232)) + (0.5869 * T_pr) + 0.1229
        elif 1.4 <= T_pr <= 2.0:
            z = (P_pr * ((0.0984 * T_pr) - 0.2053)) + (0.0621 * T_pr) + 0.8580
        elif 2.0 <= T_pr <= 3.0:
            z = (P_pr * ((0.0211 * T_pr) - 0.0527)) + (0.0127 * T_pr) + 0.9549
    elif 2.8 <= P_pr <= 5.4:
        if 1.05 <= T_pr <= 1.2:
            z = (P_pr * ((-0.3278 * T_pr) + 0.4752)) + (1.8223 * T_pr) - 1.9036
        elif 1.2 <= T_pr <= 1.4:
            z = (P_pr * ((-0.2521 * T_pr) + 0.3871)) + (1.6087 * T_pr) - 1.6635
        elif 1.4 <= T_pr <= 2.0:
            z = (P_pr * ((-0.0284 * T_pr) + 0.0625)) + (0.4714 * T_pr) - 0.0011
        elif 2.0 <= T_pr <= 3.0:
            z = (P_pr * ((0.0041 * T_pr) + 0.0039)) + (0.0607 * T_pr) + 0.7927
    elif 5.4 <= P_pr <= 15:
        if 1.05 <= T_pr <= 3.0:
            z = (P_pr * ((0.711 + (3.66 * T_pr))**(-1.4667))) - ((1.637)/((0.319 * T_pr) + 0.522)) + 2.071
    
    return z

@vectorize_decorator
def dranchuk_abou_kassem(P: float, T: float, Ppc: float, Tpc: float) -> float:
    """This is a correlation for determining the z factor for method of Dranchuk P.M., Abou-Kassem J.H.

    #### Args:
        - P (float): Pressure [psia]
        - T (float): Temperature [oR]
        - Ppc (float): Pressure pseudocritical [psia]
        - Tpc (float): Temperature pseudocritical [oR]

    #### Returns:
        z (float): Gas Compressibility Factor [dimensionless]
    """  
    from numpy import exp
    
    A1 = 0.3265
    A2 = -1.07
    A3 = -0.5339
    A4 = 0.01569
    A5 = -0.05165
    A6 = 0.5475
    A7 = -0.7361
    A8 = 0.1844
    A9 = 0.1056
    A10 = 0.6134
    A11 = 0.721
    
    P_pr = P / Ppc
    T_pr = T / Tpc
    
    z = 0.5
    tolerance = 1e-6
    max_iterations = 1000
    iteration = 0
    
    while iteration < max_iterations:
        
        dr = 0.27*(P_pr/(z*T_pr))
        
        F = z - (1 + (A1 + (A2/T_pr) + (A3/(T_pr**3)) + (A4/(T_pr**4)) + (A5/(T_pr**5)))*dr + \
        (A6 + (A7/(T_pr)) + (A8/T_pr**2))*(dr**2) - \
            A9*((A7/T_pr) + (A8/(T_pr**2)))*(dr**5) + \
                A10*(1 + A11*(dr**2))*((dr**2)/(T_pr**3))*exp(-A11*(dr**2)))
    
        dFz = 1 + (A1 + (A2/T_pr) + (A3/(T_pr**3)) + (A4/(T_pr**4)) + (A5/(T_pr**5)))*(dr/z) + \
        2*(A6 + (A7/(T_pr)) + (A8/T_pr**2))*((dr**2)/z) - \
            5*A9*((A7/T_pr) + (A8/(T_pr**2)))*((dr**5)/z) + \
                ((2*A10*(dr**2))/(z*(T_pr**3)))*(1+(A11*(dr**2))-((A11*(dr**2))**2))*exp(-A11*(dr**2))
                
        z -= F/dFz
        
        if abs(F) < tolerance:
            break
        
        iteration += 1
    
    dr = 0.27*(P_pr/(z*T_pr))
    
    z = 1 + (A1 + (A2/T_pr) + (A3/(T_pr**3)) + (A4/(T_pr**4)) + (A5/(T_pr**5)))*dr + \
        (A6 + (A7/(T_pr)) + (A8/T_pr**2))*(dr**2) - \
            A9*((A7/T_pr) + (A8/(T_pr**2)))*(dr**5) + \
                A10*(1 + A11*(dr**2))*((dr**2)/(T_pr**3))*exp(-A11*(dr**2))
                
    return z

@vectorize_decorator
def dranchuk_purvis_robinson(P: float|int, T: float|int, Ppc: float, Tpc: float) -> float:
    """This is a correlation for determining the z factor for method of Dranchuk P.M., Purvis R.A. & Robinson D.B.

    #### Args:
        - P (float): Pressure [psia]
        - T (float): Temperature [oR]
        - Ppc (float): Pressure pseudocritical [psia]
        - Tpc (float): Temperature pseudocritical [oR]

    #### Returns:
        z (float): Gas Compressibility Factor [dimensionless]
    """  
    from numpy import exp
    
    A1 = 0.31506237
    A2 = -1.0467099
    A3 = -0.57832729
    A4 = 0.53530771
    A5 = -0.61232032
    A6 = -0.10488813
    A7 = 0.68157001
    A8 = 0.68446549
    
    P_pr = P / Ppc
    T_pr = T / Tpc
    
    z = 0.5
    tolerance = 1e-6
    max_iterations = 1000
    iteration = 0
    
    while iteration < max_iterations:
        
        dr = 0.27*(P_pr/(z*T_pr))
        
        F = z - (1 + (A1 + (A2/(T_pr)) + (A3/(T_pr**3)))*dr + \
        (A4 + (A5/T_pr))*(dr**2) + ((A5*A6*(dr**5))/(T_pr)) + \
            A7*(1 + (A8*(dr**2)))*((dr**2)/(T_pr**3))*exp(-A8*(dr**2)))
        
        dFz = 1 + (A1 + (A2/T_pr) + (A3/(T_pr**3)))*(dr/z) + \
            2*(A4 + (A5/T_pr))*((dr**2)/z) + ((5*A5*A6*(dr**5))/(z*T_pr)) + \
                ((2*A7*(dr**2))/(z*(T_pr**3)))*(1 + (A8*(dr**2)) - ((A8*(dr**2))**2))*exp(-A8*(dr**2))
                
        z -= F/dFz
        
        if abs(F) < tolerance:
            break
        
        iteration += 1
    
    dr = 0.27*(P_pr/(z*T_pr))

    z = 1 + (A1 + (A2/(T_pr)) + (A3/(T_pr**3)))*dr + \
        (A4 + (A5/T_pr))*(dr**2) + ((A5*A6*(dr**5))/(T_pr)) + \
            A7*(1 + (A8*(dr**2)))*((dr**2)/(T_pr**3))*exp(-A8*(dr**2))
            
    return z


# import numpy as np
# temperatura = 194 + 460
# presion = np.arange(500, 10500, 500)
# temp_critica = 485.9
# pres_critica = 680

# print(dranchuk_purvis_robinson(presion, temperatura, pres_critica, temp_critica))
# presiones = list(range(10, 10100, 500))

# factor = [dranchuk_purvis_robinson(p, temperatura, pres_critica, temp_critica) for p in presiones]

# print(factor)

# import matplotlib.pyplot as plt

# plt.plot(presiones, factor)
# plt.show()