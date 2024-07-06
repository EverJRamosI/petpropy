""" 
## PetProPy

A comprehensive package designed to determine the physical properties of petroleum, gas, and water for the oil and gas industry. This package provides tools and functions to accurately compute essential properties such as density, viscosity, thermal conductivity, and more, helping engineers and scientists in the industry to optimize processes and make informed decisions.

### Key Features:

- Calculate physical properties of petroleum, gas, and water.
- Reliable and accurate computation methods.
- User-friendly interface for easy integration into projects.
- Detailed documentation and examples.

### Modules:

- gas: Functions to calculate the properties of natural gas.
- oil: Functions to determine the physical properties of petroleum.
- water: Functions to assess the properties of water in various conditions.
- _tools: Tools functions to support the main modules.

### Usage:

    Import the necessary modules and use the provided functions to compute physical properties. For detailed examples and usage, refer to the documentation.

### Example:

    from petpropy import z_g, mu_o
    
    z_factor = z_g.dranchuk_purvis_robinson(pressures, temperature, pres_critica, temp_critica)
    print(z_factor)
    
    muo = mu_o.mu_oil(pressures, temperature, api, rs, Pb)
    print(muo)
    
### License:

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

from petpropy.gas import z_g, c_g, rho_g, mu_g, Bg, pc_g, gamma_g, m_g
from petpropy.oil import Pb, sigma_o, c_o, rho_o, mu_o, Bo, Rs, gamma_o, Bt
from petpropy.water import sigma_w, Rsw, c_w, rho_w, mu_w, Bw
from petpropy._tools import tools

__all__ = ['z_g', 'c_g', 'rho_g', 'mu_g', 'Bg', 'pc_g', 'gamma_g', 'm_g',
            'Pb', 'sigma_o', 'c_o', 'rho_o', 'mu_o', 'Bo', 'Rs', 'gamma_o', 'Bt',
            'sigma_w', 'Rsw', 'c_w', 'rho_w', 'mu_w', 'Bw',
            'tools']