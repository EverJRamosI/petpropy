"""
### weight_molecular_gas.py
#### Functions:
    - gas_weight_molecular
    - gas_weight_molecular_df
    
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

def gas_weight_molecular(m_C7: float=0) -> float:
    """it gets of the Weight Molecular Gas. If you have the Weight Molecular C+7,\n
    you can input it into the function.  

    #### Args:
        m_C7 (float, optional): Weight Molecular of the C+7. Defaults to 0.

    #### Returns:
        m_g (float): Weight Molecular Gas [lb/lb-mol]
    """    
    m_gi = {'C1': 16.043, 'C2': 30.070, 'C3': 44.097, 'n-C4': 58.124, 'i-C4': 58.124, 'n-C5': 72.151, 'i-C5': 72.151, 'C6': 86.178, 'C7': 100.205, 'C8': 114.232, 'C9': 128.259, 'C10': 142.286, 'N2+O2': 28.963, 'CO2': 44.010, 'He': 4.003, 'H2S': 34.076, 'N2': 28.013, 'O2': 31.999, 'C7+': m_C7}
    components = list(m_gi.keys())
    
    print(f'Available mix components: {components}')
    comp_mix = [comp.strip() for comp in input('Enter the components of the mixture separated by commas: ').split(',')]
    fractions = [float(fraccion.strip()) for fraccion in input(f'Enter the fractions of the components in order {comp_mix} and separated by commas: ').split(',')]
    
    y_gi = dict(zip(comp_mix, fractions))
    
    m_gas = sum(m_gi[componente] * y_gi[componente] for componente in comp_mix)
    
    return float(m_gas)

def gas_weight_molecular_df(df_data: pd.DataFrame, m_C7: float=0) -> float:
    """it gets of the Weight Molecular Gas through the DataFrame. If you have the Weight Molecular C+7,\n
    you can input it into the function. 

    #### Args:
        df_data (pd.DataFrame): DataFrame of the mix.
        m_C7 (float, optional): Weight Molecular of the C+7. Defaults to 0.

    #### Returns:
        m_g (float): Weight Molecular Gas [lb/lb-mol]
    """    
    m_gi = {'C1': 16.043, 'C2': 30.070, 'C3': 44.097, 'n-C4': 58.124, 'i-C4': 58.124, 'n-C5': 72.151, 'i-C5': 72.151, 'C6': 86.178, 'C7': 100.205, 'C8': 114.232, 'C9': 128.259, 'C10': 142.286, 'N2+O2': 28.963, 'CO2': 44.010, 'He': 4.003, 'H2S': 34.076, 'N2': 28.013, 'O2': 31.999, 'C7+': m_C7}
    
    df_data = df_data.to_dict(orient='records')
    
    mix = [list(m.values())[0] for m in df_data]
    
    y_gi = {list(item.values())[0]: list(item.values())[1] for item in df_data}
    
    m_gas = sum(m_gi[component] * y_gi[component] for component in mix)
    
    return float(round(m_gas, 5))


