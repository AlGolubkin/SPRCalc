import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

a = 300 #0.3 micro metre #0.0000003 metres # START
b = 700 #0.7 micro metre #0.0000007 metres # STOP
step = 1 #0.001 micro metre #0.000000001 metres
N = (b - a)/step

def ref_ind_Air():
    x = np.linspace(a, b, N) # Wavelength
    n = (0.05792105 / (238.0185 - x**(-2))) + (0.00167917 / (57.362 - x**(-2))) - 1 # Dispersion formula of AIR

    n_mas = pd.DataFrame({'Wavelength' : x, 'Refractive index' : n})
    return n_mas

def ref_ind_Glass():
    x = np.linspace(a, b, N) # Wavelength
    n = (0.05792105 / (238.0185 - x**(-2))) + (0.00167917 / (57.362 - x**(-2))) - 1 # Dispersion formula of AIR

    n_mas = pd.DataFrame({'Wavelength' : x, 'Refractive index' : n})
    return n_mas

def ref_ind_Ag():
    x = np.linspace(a, b, N) # Wavelength
    n = (0.05792105 / (238.0185 - x**(-2))) + (0.00167917 / (57.362 - x**(-2))) - 1 # Dispersion formula of AIR

    n_mas = pd.DataFrame({'Wavelength' : x, 'Refractive index' : n})
    return n_mas

def ref_ind_Al2O3():
    x = np.linspace(a, b, N) # Wavelength
    n = np.sqrt(((1.4313493*x**2)/(x**2-0.0726631**2))+((0.65054713*x**2)/(x**2-0.1193242**2))+((5.3414021*x**2)/(x**2-18.028251**2))+1) # Dispersion formula of AIR

    n_mas = pd.DataFrame({'Wavelength' : x, 'Refractive index' : n})
    return n_mas

def ref_ind_Fe2O3():
    x = np.linspace(a, b, N) # Wavelength
    n = (0.05792105 / (238.0185 - x**(-2))) + (0.00167917 / (57.362 - x**(-2))) - 1 # Dispersion formula of AIR

    n_mas = pd.DataFrame({'Wavelength' : x, 'Refractive index' : n})
    return n_mas

n_mas = ref_ind_Air()

fig, ax = plt.subplots()
ax.plot(n_mas['Wavelength'], n_mas['Refractive index'], label='AIR')
plt.show()