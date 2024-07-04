"""
### pseudocritical_properties_gas.py
#### Functions:
    - mathews_roland_correlation_c7
    - sutton
    - kay
    - kay_df
    - stewart
    - stewart_df
    - brown_katz
    - brown_katz_grv
    - brown_katz_df
    
    
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

import pandas as pd

def mathews_roland_correlation_c7(m_C7: float|int, gamma_c7: float) -> dict[str, float|int]:
    """It is a correlation for determining Ppc and Tpc of the C+7 for the correlation Mathews & Roland.\n
    
    #### Args:
        - m_C7 (float | int): Weight Molecular C+7 (lbs/lb-mol)
        - gamma_c7 (float): Specific Gravity C+7

    #### Returns:
        {Ppc: value, Tpc: value} (dict): Pressure and Temperature pseudocritical [psia, oR] of the C+7. 
    """    
    from math import log10
    
    Ppc_C7 = 1188 - 431*log10((m_C7 - 61.1)) + (2319 - 852*log10((m_C7 - 53.71)))*(gamma_c7 - 0.8)
    Tpc_C7 = 608 + 364*log10((m_C7 - 71.2)) + ((2450*log10(m_C7)) - 3800)*(log10(gamma_c7))
    return {"PpcC7=": round(Ppc_C7, 5), "TpcC7=": round(Tpc_C7, 5)}

def sutton(gamma_gas: float, *, yCO2: float=0, yH2S: float=0, yN2: float=0) -> list[float]:
    """These are equations for determining Ppc and Tpc for Sutton R.P.\n
    The values for defaults in yCO2, yH2S, and yN2 are 0. However,\n
    you want to input with the values for correcting the values of the Ppc and Tpc.

    #### Args:
        - gamma_gas (float): Specific Gravity Gas [dimensionless]
        - yCO2 (float, optional): Fraction of the yCO2 [decimal]. Defaults to 0.
        - yH2S (float, optional): Fraction of the yH2S [decimal]. Defaults to 0.
        - yN2 (float, optional): Fraction of the yN2 [decimal]. Defaults to 0.

    #### Returns:
        [Ppc, Tpc] (list[float]): Pressure and Temperature pseudocritical [psia, oR]
    """    
    Ppc_gas = 756.8 - 131*gamma_gas - 3.6*(gamma_gas**2)
    Tpc_gas = 169.2 + 349.5*gamma_gas - 74*(gamma_gas**2)
    
    if yCO2 != 0 or yH2S !=0 or yN2 != 0:
        
        Ppc_gas_c = (1-yN2-yCO2-yH2S)*Ppc_gas + 493*yN2 + 1071*yCO2 + 1306*yH2S
        Tpc_gas_c = (1-yN2-yCO2-yH2S)*Tpc_gas + 227*yN2 + 548*yCO2 + 672*yH2S
        
        return [Ppc_gas_c, Tpc_gas_c]
    
    return [Ppc_gas, Tpc_gas]

def kay(Ppc7: float|int=0, Tpc7: float|int=0) -> dict[str, float]:
    """These are equations for determining Pcp and Tpc for a method of the Kay W.D.\n
    If you have the Ppc and Tpc of the C+7, you can input it into the function.
    
    #### Args:
        - Ppc7 (float, optional): Pressure Pseudocritical of the C+7. Defaults to 0.
        - Tpc7 (float, optional): Temperature Pseudocritical of the C+7. Defaults to 0.

    #### Returns:
        {Ppc: value, Tpc: value} (dict[str, float]): Pressure and Temperature pseudocritical [psia, oR]
    """
    P_pci = {'C1': 667.8, 'C2': 707.8, 'C3': 616.3, 'n-C4': 550.7, 'i-C4': 529.1, 'n-C5': 488.6, 'i-C5': 490.4, 'n-C6': 436.9, 'n-C7': 396.9, 'n-C8': 360.6, 'n-C9': 332.0, 'n-C10': 304.0, 'N2+O2': 546.9, 'CO2': 1071.0, 'He': 32.99, 'H2S': 1306.0, 'N2': 493.0, 'O2': 731.4, 'C7+': Ppc7}
    T_pci = {'C1': 343.37, 'C2': 550.09, 'C3': 666.01, 'n-C4': 765.55, 'i-C4': 734.98, 'n-C5': 845.70, 'i-C5': 829.10, 'n-C6': 913.70, 'n-C7': 972.80, 'n-C8': 1024.22, 'n-C9': 1070.68, 'n-C10': 1112.10, 'N2+O2': 238.69, 'CO2': 547.90, 'He': 9.69, 'H2S': 672.70, 'N2': 227.60, 'O2': 278.57, 'C7+': Tpc7}
    
    components = list(P_pci.keys())
    print(f'Available mix componentes: {components}\nComponentes de la mezcla disponibles: {components}')
    
    comp_mix= [comp.strip() for comp in input('Enter the components of the mixture separated by commas: \nIngrese los componentes de la mezcla separados por coma: ').split(',')]
    
    fractions = [float(fraccion.strip()) for fraccion in input(f'Enter the fractions of the components according to order {comp_mix} and separated by comma\nIntroduzca las fracciones de los componentes según el orden {comp_mix} y separado por comas: ').split(',')]
    
    y_gi = dict(zip(comp_mix, fractions))
    
    Ppc_gas = sum(P_pci[componente] * y_gi[componente] for componente in comp_mix)
    Tpc_gas = sum(T_pci[componente] * y_gi[componente] for componente in comp_mix)
    
    return {'Ppc=': Ppc_gas, 'Tpc=': Tpc_gas}

def kay_df(df_data: pd.DataFrame, Ppc7: float=0, Tpc7: float=0) -> dict[str, float]:
    """These are equations for determining Pcp and Tpc for a method of the Kay W.D.\n
    You can input a DataFrame for the function, this DataFrame have to be in the first columns the components,\n
    in the second columns the fractions of the components, including the component C+7.
    If you have the Ppc and Tpc of the C+7, you can input it into the function.
    
    #### Args:
        - df_data (pd.DataFrame): This is DataFrame with the components.
        - Ppc7 (float, optional): Pressure Pseudocritical of the C+7. Defaults to 0.
        - Tpc7 (float, optional): Temperature Pseudocritical of the C+7. Defaults to 0.

    #### Returns:
        {Ppc: value, Tpc: value} (dict[str, float]): Pressure and Temperature pseudocritical [psia, oR]
    """    
        
    P_pci = {'C1': 667.8, 'C2': 707.8, 'C3': 616.3, 'n-C4': 550.7, 'i-C4': 529.1, 'n-C5': 488.6, 'i-C5': 490.4, 'C6': 436.9, 'C7': 396.9, 'C8': 360.6, 'C9': 332.0, 'C10': 304.0, 'N2+O2': 546.9, 'CO2': 1071.0, 'He': 32.99, 'H2S': 1306.0, 'N2': 493.0, 'O2': 731.4, 'C7+': Ppc7}
    T_pci = {'C1': 343.37, 'C2': 550.09, 'C3': 666.01, 'n-C4': 765.55, 'i-C4': 734.98, 'n-C5': 845.70, 'i-C5': 829.10, 'C6': 913.70, 'C7': 972.80, 'C8': 1024.22, 'C9': 1070.68, 'C10': 1112.10, 'N2+O2': 238.69, 'CO2': 547.90, 'He': 9.69, 'H2S': 672.70, 'N2': 227.60, 'O2': 278.57, 'C7+': Tpc7}
    
    df_data = df_data.to_dict(orient='records')
    
    comp_mix = [list(m.values())[0] for m in df_data]
    
    y_gi = {list(item.values())[0]: list(item.values())[1] for item in df_data}
    
    Ppc_gas = sum(P_pci[componente] * y_gi[componente] for componente in comp_mix)
    Tpc_gas = sum(T_pci[componente] * y_gi[componente] for componente in comp_mix)
    
    return {'Ppc=': Ppc_gas, 'Tpc=': Tpc_gas}

def stewart(Ppc7: float=0, Tpc7: float=0) -> dict[str, float]:
    """These are equations for determining Pcp and Tpc for a method of the Stewart W.F. Burkhadt S.F. & Voo, D.\n
    If you have the Ppc and Tpc of the C+7, you can input it into the function.

    #### Args:
        - Ppc7 (float, optional): Pressure Pseudocritical of the C+7. Defaults to 0.
        - Tpc7 (float, optional): Temperature Pseudocritical of the C+7. Defaults to 0.

    #### Returns:
        {Ppc: value, Tpc: value} (dict[str, float]): Pressure and Temperature pseudocritical [psia, oR]
    """    
    P_pci = {'C1': 667.8, 'C2': 707.8, 'C3': 616.3, 'n-C4': 550.7, 'i-C4': 529.1, 'n-C5': 488.6, 'i-C5': 490.4, 'n-C6': 436.9, 'n-C7': 396.9, 'n-C8': 360.6, 'n-C9': 332.0, 'n-C10': 304.0, 'N2+O2': 546.9, 'CO2': 1071.0, 'He': 32.99, 'H2S': 1306.0, 'N2': 493.0, 'O2': 731.4, 'C7+': Ppc7}
    T_pci = {'C1': 343.37, 'C2': 550.09, 'C3': 666.01, 'n-C4': 765.55, 'i-C4': 734.98, 'n-C5': 845.70, 'i-C5': 829.10, 'n-C6': 913.70, 'n-C7': 972.80, 'n-C8': 1024.22, 'n-C9': 1070.68, 'n-C10': 1112.10, 'N2+O2': 238.69, 'CO2': 547.90, 'He': 9.69, 'H2S': 672.70, 'N2': 227.60, 'O2': 278.57, 'C7+': Tpc7}
    
    components = list(P_pci.keys())
    print(f'Available mix componentes: {components}\nComponentes de la mezcla disponibles: {components}')
    
    comp_mix= [comp.strip() for comp in input('Enter the components of the mixture separated by commas: \nIngrese los componentes de la mezcla separados por coma: ').split(',')]
    
    fractions = [float(fraccion.strip()) for fraccion in input(f'Enter the fractions of the components according to order {comp_mix} and separated by comma\nIntroduzca las fracciones de los componentes según el orden {comp_mix} y separado por comas: ').split(',')]
    
    y_gi = dict(zip(comp_mix, fractions))
    
    J = ((1/3) * (sum(y_gi[componente] * (T_pci[componente]/P_pci[componente]) for componente in comp_mix))) + ((2/3) * (sum(y_gi[componente] * (T_pci[componente]/P_pci[componente])**(1/2) for componente in comp_mix))**2)
    K = sum(y_gi[componente] * (T_pci[componente]/(P_pci[componente]**(1/2))) for componente in comp_mix)
    
    Tpc_gas = round((K**2)/J, 5)
    Ppc_gas = round(Tpc_gas/J, 5)
    
    return {'Ppc=': Ppc_gas, 'Tpc=': Tpc_gas}

def stewart_df(df_data: pd.DataFrame, Ppc7: float=0, Tpc7: float=0) -> dict[str, float]:
    """These are equations for determining Pcp and Tpc for a method of the Stewart W.F. Burkhadt S.F. & Voo, D for determining Ppc and Tpc.\n
    You can input a DataFrame for the function, this DataFrame have to be in the first columns the components, in the second columns the fractions of the components,\n
    including the component C+7.
    If you have the Ppc and Tpc of the C+7, you can input it into the function.

    #### Args:
        - df_data (pd.DataFrame): This is DataFrame with the components.
        - Ppc7 (float, optional): Pressure Pseudocritical of the C+7. Defaults to 0.
        - Tpc7 (float, optional): Temperature Pseudocritical of the C+7. Defaults to 0.

    #### Returns:
        {Ppc: value, Tpc: value} (dict[str, float]): Pressure and Temperature pseudocritical [psia, oR]
    """    
    P_pci = {'C1': 667.8, 'C2': 707.8, 'C3': 616.3, 'n-C4': 550.7, 'i-C4': 529.1, 'n-C5': 488.6, 'i-C5': 490.4, 'C6': 436.9, 'C7': 396.9, 'C8': 360.6, 'C9': 332.0, 'C10': 304.0, 'N2+O2': 546.9, 'CO2': 1071.0, 'He': 32.99, 'H2S': 1306.0, 'N2': 493.0, 'O2': 731.4, 'C7+': Ppc7}
    T_pci = {'C1': 343.37, 'C2': 550.09, 'C3': 666.01, 'n-C4': 765.55, 'i-C4': 734.98, 'n-C5': 845.70, 'i-C5': 829.10, 'C6': 913.70, 'C7': 972.80, 'C8': 1024.22, 'C9': 1070.68, 'C10': 1112.10, 'N2+O2': 238.69, 'CO2': 547.90, 'He': 9.69, 'H2S': 672.70, 'N2': 227.60, 'O2': 278.57, 'C7+': Tpc7}
    
    df_data = df_data.to_dict(orient='records')
    
    comp_mix = [list(m.values())[0] for m in df_data]
    
    y_gi = {list(item.values())[0]: list(item.values())[1] for item in df_data}
    
    J = ((1/3) * (sum(y_gi[componente] * (T_pci[componente]/P_pci[componente]) for componente in comp_mix))) + ((2/3) * (sum(y_gi[componente] * (T_pci[componente]/P_pci[componente])**(1/2) for componente in comp_mix))**2)
    K = sum(y_gi[componente] * (T_pci[componente]/(P_pci[componente]**(1/2))) for componente in comp_mix)
    
    Tpc_gas = round((K**2)/J, 5)
    Ppc_gas = round(Tpc_gas/J, 5)
    
    return {'Ppc=': Ppc_gas, 'Tpc=': Tpc_gas}

def brown_katz(gamma_gas: float) -> dict[str, float]:
    """These are equations for determining Pcp and Tpc for a method of the Brown G.G., Katz D.L., Oberfell G.G. & Alden R.C.\n

    #### Args:
        gamma_gas (float): Specific Gravity Gas [decimals]

    #### Returns:
        {Ppc: value, Tpc: value} (dict[str, float]): Pressure and Temperature pseudocritical [psia, oR]
    """    
    y_impurities = {'yN2': 0, 'yCO2': 0, 'yH2S': 0}
    y_impu_i = list(y_impurities.keys())
    
    f_impurity = [float(imp.strip()) for imp in input(f'Enter the following impurities {y_impu_i} in decimals in the corresponding order and separated by a comma (if you do not have impurities, enter zero [0])\nIngrese las siguiente impurezas {y_impu_i} en decimales en el orden correspondiente y separados por una coma (si no tiene una impureza coloque cero [0]): ').split(',')]
    
    y_impurities = dict(zip(y_impu_i, f_impurity))
    
    gamma_mg_hc = (gamma_gas - 0.967*y_impurities.get('yN2') - 1.52*y_impurities.get('yCO2') - 1.18*y_impurities.get('yH2S')) / (1 - y_impurities.get('yN2') - y_impurities.get('yCO2') - y_impurities.get('yH2S'))

    type_gas = input('Is the type of gas Condensed Gas? [Y/N]: \n¿El tipo de gas qué tiene es Gas Condensado? [Y/N]: ').upper()
    
    if type_gas == 'Y':  
        Ppc_hc = 706 - 51.7*gamma_mg_hc - 11.1*(gamma_mg_hc**2)
        Tpc_hc = 187 + 330*gamma_mg_hc - 71.5*(gamma_mg_hc**2)
    else:
        Ppc_hc = 677 + 15*gamma_mg_hc - 37.5*(gamma_mg_hc**2)
        Tpc_hc = 168 + 325*gamma_mg_hc - 12.5*(gamma_mg_hc**2)
        
    Ppc_gas = round((1 - y_impurities.get('yN2') - y_impurities.get('yCO2') - y_impurities.get('yH2S'))*Ppc_hc + 493*y_impurities.get('yN2') + 1071*y_impurities.get('yCO2') + 1306*y_impurities.get('yH2S'), 5)
    Tpc_gas = round((1 - y_impurities.get('yN2') - y_impurities.get('yCO2') - y_impurities.get('yH2S'))*Tpc_hc + 227*y_impurities.get('yN2') + 548*y_impurities.get('yCO2') + 672*y_impurities.get('yH2S'), 5)
    
    return {'Ppc=': Ppc_gas, 'Tpc=': Tpc_gas}

def brown_katz_grv(gamma_gas: float, *, yCO2: float=0, yH2S: float=0, yN2: float=0, type_gas: bool=True) -> list[float]:
    """These are equations for determining Pcp and Tpc for a method of the Brown G.G., Katz D.L., Oberfell G.G. & Alden R.C.\n
    The values for defaults in yCO2, yH2S, and yN2 are 0. However,\n
    you want to input with the values for correcting the values of the Ppc and Tpc.
    
    #### Args:
        - gamma_gas (float): Specific Gravity Gas [decimal]
        - yCO2 (float, optional): Fraction of the yCO2 [decimal]. Defaults to 0.
        - yH2S (float, optional): Fraction of the yH2S [decimal]. Defaults to 0.
        - yN2 (float, optional): Fraction of the yN2 [decimal]. Defaults to 0.
        - type_gas (bool, optional): The type of gas. Defaults to True for Condensed Gas and False for the Natural Gas.

    #### Returns:
        [Ppc, Tpc] (list[float]): Pressure and Temperature pseudocritical [psia, oR]
    """    
    y_ghc = (gamma_gas - 0.967*yN2 - 1.52*yCO2 - 1.18*yH2S) / (1- yN2 - yCO2 - yH2S)
    
    if type_gas == True:
        Ppc_gas = 677 + 15*y_ghc - 37.5*(y_ghc**2)
        Tpc_gas = 168 + 325*y_ghc - 12.5*(y_ghc**2)
    else:
        Ppc_gas = 706 - 51.7*y_ghc - 11.1*(y_ghc**2)
        Tpc_gas = 187 + 330*y_ghc - 71.5*(y_ghc**2)
    
    if yCO2 != 0 or yH2S !=0 or yN2 != 0:
        
        Ppc_gas_c = (1-yN2-yCO2-yH2S)*Ppc_gas + 493*yN2 + 1071*yCO2 + 1306*yH2S
        Tpc_gas_c = (1-yN2-yCO2-yH2S)*Tpc_gas + 227*yN2 + 548*yCO2 + 672*yH2S
        
        return [Ppc_gas_c, Tpc_gas_c]
    
    return [Ppc_gas, Tpc_gas]

def brown_katz_df(df_data: pd.DataFrame, gamma_gas: float) -> dict[str, float]:
    """These are equations for determining Pcp and Tpc for a method of the Brown G.G., Katz D.L., Oberfell G.G. & Alden R.C.,\n
    through DataFrame.

    #### Args:
        df_data (pd.DataFrame): This is DataFrame with the components.
        gamma_gas (float): Specific Gravity Gas [decimal]

    #### Returns:
        {Ppc: value, Tpc: value} (dict[str, float]): Pressure and Temperature pseudocritical [psia, oR]
    """    
    
    y_impurities = {'N2': 0, 'CO2': 0, 'H2S': 0}
    
    df_data = df_data.to_dict(orient='records')
        
    y_gi = {list(item.values())[0]: list(item.values())[1] for item in df_data}
        
    for key in y_impurities:
        if key in y_gi:
            y_impurities[key] = y_gi[key]
    
    gamma_mg_hc = (gamma_gas - 0.967*y_impurities.get('N2') - 1.52*y_impurities.get('CO2') - 1.18*y_impurities.get('H2S')) / (1 - y_impurities.get('N2') - y_impurities.get('CO2') - y_impurities.get('H2S'))
    
    type_gas = input('Is the type of gas Condensed Gas? [Y/N]: \n¿El tipo de gas qué tiene es Gas Condensado? [Y/N]: ').upper()
    
    if type_gas == 'Y':  
        Ppc_hc = 706 - 51.7*gamma_mg_hc - 11.1*(gamma_mg_hc**2)
        Tpc_hc = 187 + 330*gamma_mg_hc - 71.5*(gamma_mg_hc**2)
    else:
        Ppc_hc = 677 + 15*gamma_mg_hc - 37.5*(gamma_mg_hc**2)
        Tpc_hc = 168 + 325*gamma_mg_hc - 12.5*(gamma_mg_hc**2)
    
    Ppc_gas = round((1 - y_impurities.get('N2') - y_impurities.get('CO2') - y_impurities.get('H2S'))*Ppc_hc + 493*y_impurities.get('N2') + 1071*y_impurities.get('CO2') + 1306*y_impurities.get('H2S'), 5)
    Tpc_gas = round((1 - y_impurities.get('N2') - y_impurities.get('CO2') - y_impurities.get('H2S'))*Tpc_hc + 227*y_impurities.get('N2') + 548*y_impurities.get('CO2') + 672*y_impurities.get('H2S'), 5)
        
    return {'Ppc=': Ppc_gas, 'Tpc=': Tpc_gas}

# if __name__ == '__main__':
#     pass


