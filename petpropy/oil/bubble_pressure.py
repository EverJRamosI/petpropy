"""
### bubble_pressure.py
#### Functions:
    - standing
    - lasater
    - vazquez_beggs
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
@email: everramosisla@gmail.com
"""

def standing(R_sb: float, gamma_gas: float, T: float|int, gamma_api: float|int, *, y_N2: float=0, y_CO2: float=0, y_H2S: float=0) -> float:
    """This is a correlation for determining the Bubble Point Pressure, through Standing M.B.\n
    If you want to correct the bubble pressure, it inputs the fractions of the components.
    
    #### Args:
        R_sb (float): Solution Gas-Oil Ratio to P>=Pb [PCN/BN]
        gamma_gas (float): Specific Gravity Gas [fraction]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        y_N2 (float, optional): Fraction of the yN2 [decimal]. Defaults to 0.
        y_CO2 (float, optional): Fraction of the yCO2 [decimal]. Defaults to 0.
        y_H2S (float, optional): Fraction of the yH2S [decimal]. Defaults to 0.

    #### Returns:
        Pb (float): Bubble Point Pressure [psia]
    """    
    F = (R_sb / gamma_gas)**0.83 * (10**(0.00091*(T-460) - 0.0125*gamma_api))
    Pb = 18.2 * (F -1.4)
    
    if y_N2 != 0 or y_CO2 != 0 or y_H2S != 0:
        
        T -= 460
        c_N2 = 1 + (((((-2.65e-4*gamma_api)+(5.5e-3))*T)+((0.0931*gamma_api)-(0.8295)))*y_N2) + ((((1.954e-11*(gamma_api**4.699))*T)+((0.027*gamma_api)-(2.366)))*(y_N2**2))
        c_CO2 = 1 - (693.8*y_CO2*(T**(-1.553)))
        c_H2S = 1 - ((0.9035+(0.0015*gamma_api))*y_H2S) + ((0.019*(45-gamma_api))*(y_H2S**2))
        Pb_c = Pb*c_N2*c_CO2*c_H2S
        
        return round(Pb_c, 5)
        
    return round(Pb, 5)

def lasater(R_sb: float, gamma_gas: float, T: float|int, gamma_api: float|int, *, y_N2: float=0, y_CO2: float=0, y_H2S: float=0) -> float:
    """This is a correlation for determining the Bubble Point Pressure, through Lasater J.A.\n
    If you want to correct the bubble pressure, it inputs the fractions of the components. 
    
    #### Args:
        R_sb (float): Solution Gas-Oil Ratio to P>=Pb [PCN/BN]
        gamma_gas (float): Specific Gravity Gas [fraction]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        y_N2 (float, optional): Fraction of the yN2 [decimal]. Defaults to 0.
        y_CO2 (float, optional): Fraction of the yCO2 [decimal]. Defaults to 0.
        y_H2S (float, optional): Fraction of the yH2S [decimal]. Defaults to 0.

    #### Returns:
        Pb (float): Bubble Point Pressure [psia]
    """    
    from numpy import exp
    
    if gamma_api <= 40:
        m_oil = 630 - (10 * gamma_api)
    elif gamma_api > 40:
        m_oil = 73110 / (gamma_api**1.562)
    
    gamma_oil = 141.5 / (gamma_api + 131.5)
    
    y_gas = (R_sb/379.3) / ((R_sb/379.3) + ((350*gamma_oil)/(m_oil)))
    
    if y_gas <= 0.60:
        Pf = (0.679*exp(2.786*y_gas)) - 0.323
    elif y_gas > 0.60:
        Pf = (8.26*(y_gas**3.56)) + 1.95
    
    Pb = Pf * (T/gamma_gas)
    
    if y_N2 != 0 or y_CO2 != 0 or y_H2S != 0:
        
        T -= 460
        c_N2 = 1 + (((((-2.65e-4*gamma_api)+(5.5e-3))*T)+((0.0931*gamma_api)-(0.8295)))*y_N2) + ((((1.954e-11*(gamma_api**4.699))*T)+((0.027*gamma_api)-(2.366)))*(y_N2**2))
        c_CO2 = 1 - (693.8*y_CO2*(T**(-1.553)))
        c_H2S = 1 - ((0.9035+(0.0015*gamma_api))*y_H2S) + ((0.019*(45-gamma_api))*(y_H2S**2))
        Pb_c = Pb*c_N2*c_CO2*c_H2S
        
        return round(Pb_c, 5)
    
    return round(Pb, 5)

def vazquez_beggs(R_sb: float, gamma_gas: float, T: float|int, gamma_api: float|int, *, P_sp: float=0, T_sp: float=0, y_N2: float=0, y_CO2: float=0, y_H2S: float=0) -> float:
    """This is a correlation for determining the Bubble Point Pressure, through Vazquez M.E. & Beggs H.D.\n
    If you want to correct the bubble pressure, it inputs the fractions of the components. 
    
    #### Args:
        R_sb (float): Solution Gas-Oil Ratio to P>=Pb [PCN/BN]
        gamma_gas (float): Specific Gravity Gas [fraction]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        P_sp (float, optional): Separator Pressure [psia]. Defaults to 0.
        T_sp (float, optional): Separator Temperature [oR]. Defaults to 0.
        y_N2 (float, optional): Fraction of the yN2 [decimal]. Defaults to 0.
        y_CO2 (float, optional): Fraction of the yCO2 [decimal]. Defaults to 0.
        y_H2S (float, optional): Fraction of the yH2S [decimal]. Defaults to 0.

    #### Returns:
        Pb (float): Bubble Point Pressure [psia]
    """    
    from numpy import exp, log10
        
    if gamma_api <= 30:
        C1 = 0.0362
        C2 = 1.0937
        C3 = 25.724
    elif gamma_api > 30:
        C1 = 0.0178
        C2 = 1.1870
        C3 = 23.931
        
    if P_sp != 0 and T_sp != 0:
        T -= 460
        T_sp -= 460
        gamma_gc = gamma_gas * (1 + (5.912e-5 * gamma_api * T_sp * log10(P_sp/114.7)))
        Pb = (R_sb/(C1 * gamma_gc * exp((C3 * gamma_api)/(T))))**(1/C2)
    else:
        T -= 460
        Pb = (R_sb/(C1 * gamma_gas * exp((C3 * gamma_api)/(T))))**(1/C2)
    
    if y_N2 != 0 or y_CO2 != 0 or y_H2S != 0:
        T -= 460
        c_N2 = 1 + (((((-2.65e-4*gamma_api)+(5.5e-3))*T)+((0.0931*gamma_api)-(0.8295)))*y_N2) + ((((1.954e-11*(gamma_api**4.699))*T)+((0.027*gamma_api)-(2.366)))*(y_N2**2))
        c_CO2 = 1 - (693.8*y_CO2*(T**(-1.553)))
        c_H2S = 1 - ((0.9035+(0.0015*gamma_api))*y_H2S) + ((0.019*(45-gamma_api))*(y_H2S**2))
        Pb_c = Pb*c_N2*c_CO2*c_H2S
        
        return round(Pb_c, 5)
    
    return round(Pb, 5)

def glaso(R_sb: float, gamma_gas: float, T: float|int, gamma_api: float|int, *, y_N2: float=0, y_CO2: float=0, y_H2S: float=0) -> float:
    """This is a correlation for determining the Bubble Point Pressure, through Glaso O.\n
    If you want to correct the bubble pressure, it inputs the fractions of the components.
    
    #### Args:
        R_sb (float): Solution Gas-Oil Ratio to P>=Pb [PCN/BN]
        gamma_gas (float): Specific Gravity Gas [fraction]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        y_N2 (float, optional): Fraction of the yN2 [decimal]. Defaults to 0.
        y_CO2 (float, optional): Fraction of the yCO2 [decimal]. Defaults to 0.
        y_H2S (float, optional): Fraction of the yH2S [decimal]. Defaults to 0.

    #### Returns:
        Pb (float): Bubble Point Pressure [psia]
    """ 
    from numpy import log10
    
    F = ((R_sb/gamma_gas)**0.816) * (((T-460)**0.172) / (gamma_api**0.989))
    
    Pb = 10**(1.7669 + (1.7447 * log10(F)) - (0.30218 * ((log10(F))**2)))
    
    if y_N2 != 0 or y_CO2 != 0 or y_H2S != 0:
        
        T -= 460
        c_N2 = 1 + (((((-2.65e-4*gamma_api)+(5.5e-3))*T)+((0.0931*gamma_api)-(0.8295)))*y_N2) + ((((1.954e-11*(gamma_api**4.699))*T)+((0.027*gamma_api)-(2.366)))*(y_N2**2))
        c_CO2 = 1 - (693.8*y_CO2*(T**(-1.553)))
        c_H2S = 1 - ((0.9035+(0.0015*gamma_api))*y_H2S) + ((0.019*(45-gamma_api))*(y_H2S**2))
        Pb_c = Pb*c_N2*c_CO2*c_H2S
        
        return round(Pb_c, 5)
    
    return round(Pb, 5)

def total(R_sb: float, gamma_gas: float, T: float|int, gamma_api: float|int, *, y_N2: float=0, y_CO2: float=0, y_H2S: float=0) -> float:
    """This is a correlation for determining the Bubble Point Pressure, through TOTAL C.F.P. J.A.\n
    If you want to correct the bubble pressure, it inputs the fractions of the components.
    
    #### Args:
        R_sb (float): Solution Gas-Oil Ratio to P>=Pb [PCN/BN]
        gamma_gas (float): Specific Gravity Gas [fraction]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        y_N2 (float, optional): Fraction of the yN2 [decimal]. Defaults to 0.
        y_CO2 (float, optional): Fraction of the yCO2 [decimal]. Defaults to 0.
        y_H2S (float, optional): Fraction of the yH2S [decimal]. Defaults to 0.

    #### Returns:
        Pb (float): Bubble Point Pressure [psia]
    """ 
    if gamma_api <= 10:
        C1 = 12.847
        C2 = 0.9636
        C3 = 0.000993
        C4 = 0.034170
    elif 10 < gamma_api <= 35:
        C1 = 25.2755
        C2 = 0.7617
        C3 = 0.000835
        C4 = 0.011292
    elif 35 < gamma_api <= 45:
        C1 = 216.4711
        C2 = 0.6922
        C3 = -0.000427
        C4 = 0.023140
    
    Pb = C1 * ((R_sb/gamma_gas)**C2) * (10**(C3*(T-460) - C4*gamma_api))
    
    if y_N2 != 0 or y_CO2 != 0 or y_H2S != 0:
        
        T -= 460
        c_N2 = 1 + (((((-2.65e-4*gamma_api)+(5.5e-3))*T)+((0.0931*gamma_api)-(0.8295)))*y_N2) + ((((1.954e-11*(gamma_api**4.699))*T)+((0.027*gamma_api)-(2.366)))*(y_N2**2))
        c_CO2 = 1 - (693.8*y_CO2*(T**(-1.553)))
        c_H2S = 1 - ((0.9035+(0.0015*gamma_api))*y_H2S) + ((0.019*(45-gamma_api))*(y_H2S**2))
        Pb_c = Pb*c_N2*c_CO2*c_H2S
        
        return round(Pb_c, 5)
    
    return round(Pb, 5)

def almarhoun(R_sb: float, gamma_gas: float, T: float|int, gamma_api: float|int, *, y_N2: float=0, y_CO2: float=0, y_H2S: float=0) -> float:
    """This is a correlation for determining the Bubble Point Pressure, through Al-Marhoun M.A.\n
    If you want to correct the bubble pressure, it inputs the fractions of the components.
    
    #### Args:
        R_sb (float): Solution Gas-Oil Ratio to P>=Pb [PCN/BN]
        gamma_gas (float): Specific Gravity Gas [fraction]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        y_N2 (float, optional): Fraction of the yN2 [decimal]. Defaults to 0.
        y_CO2 (float, optional): Fraction of the yCO2 [decimal]. Defaults to 0.
        y_H2S (float, optional): Fraction of the yH2S [decimal]. Defaults to 0.

    #### Returns:
        Pb (float): Bubble Point Pressure [psia]
    """ 
    gamma_oil = 141.5 / (gamma_api + 131.5)
    
    Pb = 5.38088e-3 * (R_sb**0.715082) * (gamma_gas**(-1.87784)) * (gamma_oil**3.1437) * (T**1.32657)
    
    if y_N2 != 0 or y_CO2 != 0 or y_H2S != 0:
        
        T -= 460
        c_N2 = 1 + (((((-2.65e-4*gamma_api)+(5.5e-3))*T)+((0.0931*gamma_api)-(0.8295)))*y_N2) + ((((1.954e-11*(gamma_api**4.699))*T)+((0.027*gamma_api)-(2.366)))*(y_N2**2))
        c_CO2 = 1 - (693.8*y_CO2*(T**(-1.553)))
        c_H2S = 1 - ((0.9035+(0.0015*gamma_api))*y_H2S) + ((0.019*(45-gamma_api))*(y_H2S**2))
        Pb_c = Pb*c_N2*c_CO2*c_H2S
        
        return round(Pb_c, 5)
    
    return round(Pb, 5)

def dokla_osman(R_sb: float, gamma_gas: float, T: float|int, gamma_api: float|int, *, y_N2: float=0, y_CO2: float=0, y_H2S: float=0) -> float:
    """This is a correlation for determining the Bubble Point Pressure, through Dokla M.E. & Osman M.E.\n
    If you want to correct the bubble pressure, it inputs the fractions of the components. 
    
    #### Args:
        R_sb (float): Solution Gas-Oil Ratio to P>=Pb [PCN/BN]
        gamma_gas (float): Specific Gravity Gas [fraction]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        y_N2 (float, optional): Fraction of the yN2 [decimal]. Defaults to 0.
        y_CO2 (float, optional): Fraction of the yCO2 [decimal]. Defaults to 0.
        y_H2S (float, optional): Fraction of the yH2S [decimal]. Defaults to 0.

    #### Returns:
        Pb (float): Bubble Point Pressure [psia]
    """ 
    gamma_oil = 141.5 / (gamma_api + 131.5)
    
    Pb = 0.836386e4 * (R_sb**0.724047) * (gamma_gas**(-1.01049)) * (gamma_oil**0.107991) * (T**(-0.952584))
    
    if y_N2 != 0 or y_CO2 != 0 or y_H2S != 0:
        
        T -= 460
        c_N2 = 1 + (((((-2.65e-4*gamma_api)+(5.5e-3))*T)+((0.0931*gamma_api)-(0.8295)))*y_N2) + ((((1.954e-11*(gamma_api**4.699))*T)+((0.027*gamma_api)-(2.366)))*(y_N2**2))
        c_CO2 = 1 - (693.8*y_CO2*(T**(-1.553)))
        c_H2S = 1 - ((0.9035+(0.0015*gamma_api))*y_H2S) + ((0.019*(45-gamma_api))*(y_H2S**2))
        Pb_c = Pb*c_N2*c_CO2*c_H2S
        
        return round(Pb_c, 5)
    
    return round(Pb, 5)

def petrosky_farshad(R_sb: float, gamma_gas: float, T: float|int, gamma_api: float|int, *, y_N2: float=0, y_CO2: float=0, y_H2S: float=0) -> float:
    """This is a correlation for determining the Bubble Point Pressure, through Petrosky G.E. & Farshad F.F.\n
    If you want to correct the bubble pressure, it inputs the fractions of the components. 
    
    #### Args:
        R_sb (float): Solution Gas-Oil Ratio to P>=Pb [PCN/BN]
        gamma_gas (float): Specific Gravity Gas [fraction]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        y_N2 (float, optional): Fraction of the yN2 [decimal]. Defaults to 0.
        y_CO2 (float, optional): Fraction of the yCO2 [decimal]. Defaults to 0.
        y_H2S (float, optional): Fraction of the yH2S [decimal]. Defaults to 0.

    #### Returns:
        Pb (float): Bubble Point Pressure [psia]
    """ 
    F = ((R_sb**0.5774)/(gamma_gas**0.8439)) * (10**((4.561e-5*((T-460)**1.3911)) - (7.916e-4*(gamma_api**1.541))))
    
    Pb = 112.727 * (F - 12.34)
    
    if y_N2 != 0 or y_CO2 != 0 or y_H2S != 0:
        
        T -= 460
        c_N2 = 1 + (((((-2.65e-4*gamma_api)+(5.5e-3))*T)+((0.0931*gamma_api)-(0.8295)))*y_N2) + ((((1.954e-11*(gamma_api**4.699))*T)+((0.027*gamma_api)-(2.366)))*(y_N2**2))
        c_CO2 = 1 - (693.8*y_CO2*(T**(-1.553)))
        c_H2S = 1 - ((0.9035+(0.0015*gamma_api))*y_H2S) + ((0.019*(45-gamma_api))*(y_H2S**2))
        Pb_c = Pb*c_N2*c_CO2*c_H2S
        
        return round(Pb_c, 5)
    
    return round(Pb, 5)

def kartoatmodjo_schmidt(R_sb: float, gamma_gas: float, T: float|int, gamma_api: float|int, *, P_sp: float=0, T_sp: float=0, y_N2: float=0, y_CO2: float=0, y_H2S: float=0):
    """This is a correlation for determining the Bubble Point Pressure, through Kartoatmodjo T. & Schimidt Z.\n
    If you want to correct the bubble pressure, it inputs the fractions of the components. 
    
    #### Args:
        R_sb (float): Solution Gas-Oil Ratio to P>=Pb [PCN/BN]
        gamma_gas (float): Specific Gravity Gas [fraction]
        T (float | int): Temperature [oR]
        gamma_api (float | int): API Gravity
        P_sp (float, optional): Separator Pressure [psia]. Defaults to 0.
        T_sp (float, optional): Separator Temperature [oR]. Defaults to 0.
        y_N2 (float, optional): Fraction of the yN2 [decimal]. Defaults to 0.
        y_CO2 (float, optional): Fraction of the yCO2 [decimal]. Defaults to 0.
        y_H2S (float, optional): Fraction of the yH2S [decimal]. Defaults to 0.

    #### Returns:
        Pb (float): Bubble Point Pressure [psia]
    """  
    from numpy import log10
            
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
    
    if P_sp != 0 and T_sp != 0:
        T -= 460
        T_sp -= 460
        gamma_gc = gamma_gas * (1 + (0.1595 * (gamma_api**0.4078) * (T_sp**(-0.2466)) * log10(P_sp/114.7))) 
        Pb = (R_sb/(C1 * (gamma_gc**C2) * (10**((C3*gamma_api)/(T)))))**C4
    else:
        T -= 460
        Pb = (R_sb/(C1 * (gamma_gas**C2) * (10**((C3*gamma_api)/(T)))))**C4
    
    if y_N2 != 0 or y_CO2 != 0 or y_H2S != 0:
        T -= 460
        c_N2 = 1 + (((((-2.65e-4*gamma_api)+(5.5e-3))*T)+((0.0931*gamma_api)-(0.8295)))*y_N2) + ((((1.954e-11*(gamma_api**4.699))*T)+((0.027*gamma_api)-(2.366)))*(y_N2**2))
        c_CO2 = 1 - (693.8*y_CO2*(T**(-1.553)))
        c_H2S = 1 - ((0.9035+(0.0015*gamma_api))*y_H2S) + ((0.019*(45-gamma_api))*(y_H2S**2))
        Pb_c = Pb*c_N2*c_CO2*c_H2S
        
        return round(Pb_c, 5)
    
    return round(Pb, 5)

# print("Standing", standing(675, 0.95, (460+180), 31, y_CO2=0.2, y_H2S=0.1))
# print("Lastar", lasater(675, 0.95, (460+180), 31, y_CO2=0.2, y_H2S=0.1))
# print("Vazquez",vazquez_beggs(675, 0.95, (460+180), 31, P_sp=100, T_sp=85, y_CO2=0.2, y_H2S=0.1))
# print("Glaso",glaso(675, 0.95, (460+180), 31, y_CO2=0.2, y_H2S=0.1))
# print("Total",total(675, 0.95, (460+180), 31, y_CO2=0.2, y_H2S=0.1))
# print("Al Marhoun",almarhoun(675, 0.95, (460+180), 31, y_CO2=0.2, y_H2S=0.1))
# print("Dokla",dokla_osman(675, 0.95, (460+180), 31, y_CO2=0.2, y_H2S=0.1))
# print("Petrosky",petrosky_farshad(675, 0.95, (460+180), 31, y_CO2=0.2, y_H2S=0.1))
# print("Kartoatmok",kartoatmodjo_schmidt(675, 0.95, (460+180), 31, P_sp=100, T_sp=85, y_CO2=0.2, y_H2S=0.1))
