"""
### solution_gas_water_ratio.py
#### Functions:
    - culberson_macketta
    - mccoy
    
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
def culberson_macketta(P: float|int, T: float|int, *, S: int=10000) -> float:
    """This is a correlation for determining the solution gas-water ratio, through Culberson O.L. & McKetta J.J.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        S (int, optional): Salinity [ppm]. Defaults to 10000.

    #### Returns:
        Rsw (float): Solution Gas-Water Ratio [PCN/STB]
    """    
    S /= 10000
    
    A = 8.15839 - (6.12265e-2*(T-460)) + (1.91663e-4*((T-460)**2)) - (2.1654e-7*((T-460)**3))
    B = (1.01021e-2) - (7.44121e-5*(T-460)) + (3.05553e-7*((T-460)**2)) - (2.94883e-10*((T-460)**3))
    C = ((-9.02505) + (0.130237*(T-460)) - (8.53425e-4*((T-460)**2)) + (2.34122e-6*((T-460)**3)) - (2.37049e-9*((T-460)**4))) * 1e-7
    
    R_swp = A + (B*P) + (C*(P**2))
    
    Rsw = (10**((-0.0840655)*S*((T-460)**(-0.285854)))) * R_swp
    
    return round(Rsw, 4)

@vectorize_decorator
def mccoy(P: float|int, T: float|int, *, S: int=10000) -> float:
    """This is a correlation for determining the solution gas-water ratio, through McCoy R.L.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        S (int, optional): Salinity [ppm]. Defaults to 10000.

    #### Returns:
        Rsw (float): Solution Gas-Water Ratio [PCN/STB]
    """    
    S /= 10000
    
    A = 2.12 + (3.45e-3*(T-460)) - (3.59e-5*((T-460)**2))
    B = 0.0107 - (5.26e-5*(T-460)) + (1.48e-7*((T-460)**2))
    C = (-8.75e-7) + (3.9e-9*(T-460)) - (1.02e-11*((T-460)**2))
    
    R_swp = A + (B*P) + (C*(P**2))
    
    Rsw = (1 - ((0.0753-(1.73e-4*(T-460)))*S)) * R_swp
    
    return round(Rsw, 4)

#print('culberson', culberson_macketta(5000, (200+460), S=20000))
#print('mccoy', mccoy(5000, (200+460), S=20000))

# import matplotlib.pyplot as plt
# presiones = list(range(500, 5100, 100))
# temperatura = 200 + 460
# #solubilidad = [culberson_macketta(p, temperatura, S=20000) for p in presiones]
# solubilidad = [mccoy(p, temperatura, S=20000) for p in presiones]
# print(presiones)
# print(solubilidad)
# plt.plot(presiones, solubilidad)
# plt.show()