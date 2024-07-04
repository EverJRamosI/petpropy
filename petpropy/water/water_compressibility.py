"""
### water_compressibility.py
#### Functions:
    - dodson_standing
    - osif
    - brill_beggs
    
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
def dodson_standing(P: float|int, T: float|int, Rsw: float, *, S: int=10000) -> float:
    """This is a correlation for determining the water compressibility, through Dodson C.R. & Standing M.B.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        Rsw (float): Solution Gas-Water Ratio [SCF/STB]
        S (int, optional): Salinity [ppm]. Defaults to 10000.

    #### Returns:
        cw (float): Water Compressibility [psia^-1]
    """    
    T -= 460
    S /= 10000
    
    A = 3.8546 - (1.34e-4*P)
    B = (-0.01052) + (4.77e-7*P)
    C = (3.9267e-5) - (8.8e-10*P)
    
    c_wp = (A + (B*T) + (C*(T**2))) / (10**6)
    
    R1c_wp = 1 + (8.9e-3*Rsw)
    
    cw1 = c_wp*R1c_wp
    
    R2c_wp = (1 + ((S**0.7)*((-5.2e-2)+(2.7e-4*T)-(1.14e-6*(T**2))+(1.121e-9*(T**3)))))
    
    cw = cw1*R2c_wp
    
    return cw

@vectorize_decorator
def osif(P: float|int, T: float|int, *, S: int=10000) -> float:
    """This is a correlation for determining the water compressibility, through Osif T.L.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        S (int, optional): Salinity [ppm]. Defaults to 10000.

    #### Returns:
        cw (float): Water Compressibility [psia^-1]
    """   
    T -= 460
    S /= 10000
    
    cw = 1 / ((7.033*P) + (541.5*S) - (537*T) + 403300)
    
    return cw

@vectorize_decorator
def brill_beggs(P: float|int, T: float|int) -> float:
    """This is a correlation for determining the water compressibility, through Brill & Beggs.

    #### Args:
        P (float | int): Pressure [psia]
        T (float | int): Temperature [oR]
        
    #### Returns:
        cw (float): Water Compressibility [psia^-1]
    """ 
    T -= 460
    
    C1 = 3.8546 - (0.000134*P)
    C2 = (-0.01052) + (4.77e-7*P)
    C3 = 3.9267e-5 - (8.8e-10*P)
    
    cw = (C1 + (C2*T) + (C3*(T**2)))*1e-6
    
    return cw
    
#print('dodson', dodson_standing(5000, (200+460), 17.8, S=20000))
#print('osif', osif(5000, (200+460), S=20000))

#import matplotlib.pyplot as plt
#from solution_gas_water_ratio import mccoy

#presiones = list(range(500, 5100, 100))
#
#temperatura = 200 + 460
#
#disolubilidad = mccoy(presiones, temperatura, S=20000)
#
#compresibilidad = osif(presiones, temperatura, S=20000)
# compresibilidad = [dodson_standing(p, temperatura, rs, S=20000) for p, rs in zip(presiones, disolubilidad)]
# #compresibilidad = [brill_beggs(p, temperatura) for p in presiones]

#print(presiones)
#print(disolubilidad)
#print(compresibilidad)
#
#plt.plot(presiones, compresibilidad)
#plt.show()