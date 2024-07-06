<h1 align="center">
<img src="https://raw.githubusercontent.com/EverJRamosI/petpropy/main/docs/images/PETPROPY.png" width="250">
</h1>

![stars](https://img.shields.io/github/stars/EverJRamosI/petpropy)
![forks](https://img.shields.io/github/forks/EverJRamosI/petpropy)


# PetProPy 

PetProPy is the fundamental package for petroleum, chemical engineering, and other similar industries.
This package is useful for determining the Physical Properties of Gas, Oil, and Water.

## Table of Contents ##

- [Installation](#installation)
- [Usage](#usage)
- [Theory](#theory)
- [Dependencies](#dependencies)
- [License](#license)

## Installation
You can install [PetProPy](https://github.com/EverJRamosI/petpropy/tree/main/petpropy) using pip:
```bash
pip install petpropy
```

## Usage
Here is a quick example of how to use PetProPy to calculate the compressibility factor of gas at a given temperature and pressure:

**INPUT**
```python
import matplotlib.pyplot as plt
import numpy as np
from petpropy import z_g

pressures = np.arange(500, 10500, 500)
temperature = 654
temp_critica = 485.9
pres_critica = 680

z_factor = z_g.dranchuk_purvis_robinson(pressures, temperature,pres_critica, temp_critica)
print(z_factor)

plt.plot(pressures, z_factor)
plt.xlabel('Pressures')
plt.ylabel('Z Factor')
plt.show()
```
**OUTPUT**
```bash
[0.89438476 0.78880765 0.7035361  0.6662592  0.67374447 0.70713281
 0.75391071 0.80757742 0.86484302 0.92400223 0.98413323 1.04471578
 1.10544521 1.16613807 1.22668188 1.28700724 1.34707163 1.40684971
 1.46632734 1.52549774]
```
![Z Factor](https://raw.githubusercontent.com/EverJRamosI/petpropy/main/docs/images/output.png)
> Pressures vs Z Factor (***[FOR MORE INFORMATION](book_test.ipynb)***)

## Theory
In the analysis of reservoir behavior, calculation of reserves and equipment design, knowledge of the physical properties of fluids is required.
The set of tests necessary to determine these properties is called Pressure-Volume-Temperature (PVT) analysis, and consists of determining the relationships between pressure, volume and temperature for a particular mixture of hydrocarbons (liquids and gases).

<details>
<summary>
Click to expand the theory section
</summary>

### Gas Properties
The properties of gases are calculated using the different correlations, for example:
- **Weight Molecular Gas**
- **Specific Gravity Gas**
- **Pseudocritical Properties Gas**
  - Kay W.D.
  - Stewart W.F., Burkhardt S.F. & Voo D.
  - Brown G.G., Katz D.L., Oberfell G.G. & Alden R.C.
- **Gas Compressibility Factor Z**
  - Correction of Whichert & Aziz
  - Standing M.B. & Katz D.L.
  - Papay J.
  - Hall K.R. & Yarborough L.
  - Brill J.P. & Beggs H.D.
  - Dranchuk P.M., Purvis R.A. & Robinson D.B
  - Dranchuk P.M. & Abou-Kassem J.H.
  - Gopal V.N.
- **Gas Volumetric Factor**
- **Gas Compressibility**
- **Gas Viscosity**
  - Lee A.L., Gonzalez M.H. & Eakin B.E.
- **Gas Density**  
### Oil Properties
The properties of oils are calculated using the different correlations, for example:
- **Specific Gravity Oil**
- **Bubble Pressure**
  - Standing M.B.
  - Lasater J.A.
  - Vazquez M.E. & Beggs H.D.
  - Glaso O.
  - TOTAL C.F.P.
  - Al-Marhoun M.A.
  - Dokla M.E. & Osman M.E.
  - Petrosky G.E. & Farshad F.F.
  - Kartoatmodjo T. & Schmidt Z.
- **Solution Gas-Oil Ratio**
  - Standing M.B.
  - Lasater J.A.
  - Vazquez M.E. & Beggs H.D.
  - Glaso O.
  - TOTAL C.F.P.
  - Al-Marhoun M.A.
  - Dokla Petrosky G.E. & Farshad F.F.
  - Kartoatmodjo T. & Schmidt Z.
- **Oil Volumetric Factor**
  - Standing M.B.
  - Vazquez M.E. & Beggs H.D.
  - Glaso O.
  - TOTAL C.F.P.
  - Al-Marhoun M.A.
  - Dokla M.E. & Osman M.E.
  - Petrosky G.E. & Farshad F.F.
  - Kartoatmodjo T. & Schmidt Z.
- **Total Volumetric Factor**
  - Glaso O.
  - Al-Marhoun M.A.
- **Oil Compressibility**
  - Vazquez M.E. & Beggs H.D.
  - Petrosky G.E. & Farshad F.F.
  - Kartoatmodjo T. & Schmidt Z.
  - McCain W.D., Rollins J.B. & Villena-Lanzi A.J.
- **Oil Viscosity**
  - **Dead Oil Viscosity**
    - Beal C.
    - Beggs H.B. & Robinson J.R.
    - Glaso O.
    - Egbogah E.O.
    - Kartoatmodjo T. & Schmidt Z.
  - **Saturated Oil Viscosity**
    - Chew J.N. & Connally C.A.
    - Beggs H.D. & Robinson J.R.
    - Kartoatmodjo T. & Schmidt Z.
  - **Unsaturated Oil Viscosity**
    - Beal C.
    - Vazquez M.E. & Beggs H.D.
    - Kartoatmodjo T. & Schmidt Z.
- **Oil Density**
- **Gas-Oil Interfacial Tension**
  - Baker O. & Swerdloff W.
### Water Properties
- **Solution Gas-Water Ratio**
  - Culberson O.L. & McKetta J.J.
  - McCoy R.L.
- **Water Volumetric Factor**
  - McCain W.D.
  - McCoy R.L.
- **Water Compressibility**
  - Dodson C.R. & Standing M.B.
  - Osif T.L.
  - Ramey H.J.
- **Water Viscosity**
  - Van Wingen N.
  - Matthews C.S. & Russel D.G.
  - McCain W.D.
  - McCoy R.L.
- **Water Density**
  - McCain W.D.
- **Gas-Water Interfacial Tension**
  - Jennings H.Y. & Newman G.H.

</details>

## Dependencies
PetProPy requires the following Python (=>3.12) libraries:
1. NumPy
2. Pandas 

```bash
pip install -r requirements.txt
```

## License
This project is licensed under the GNU Lesser General Public License v3.0 - see the [LICENSE](https://github.com/EverJRamosI/petpropy/blob/main/LICENSE) file for details.

## API Documentation
For detailed API documentation, please refer to the [official documentation](https://github.com/EverJRamosI/petpropy).