import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

class attenuation_fac:
    def __init__(self, cross_area):
        # cross_area here is a function 
        self.cross_area = cross_area
        self.h = 6.626E-34
    
    def energy_to_wavelength(self, energy, mass):
        if energy == 0: 
            wavelength = np.sqrt(self.h **2/(2 *mass * (energy+0.0000001)))
        else:           
            wavelength = np.sqrt(self.h **2/(2 *mass * (energy)))
        return wavelength
    
    # interpolate the neutron cross sections for elements
    def cross_section_interpolate(self, cross_section_sigma, cross_section_energy, mass):
        cross_section_wavelength = []
        for i in range(len(cross_section_energy)):
            wavelength = self.energy_to_wavelength(cross_section_energy[i], mass)
            cross_section_wavelength.append(wavelength)

        # sort wavelength and sigma
        cross_section_wavelength, cross_section_sigma =zip(*sorted(zip(cross_section_wavelength, cross_section_sigma)))
        print('ross', cross_section_wavelength)
        print(cross_section_sigma)
        spl = CubicSpline(cross_section_wavelength, cross_section_sigma)
        return spl

    def factor(self, molmass, mass , distance, cross_section_sigma, cross_section_energy):
        n = mass/molmass
        cross_section_wavelength = []
        for i in range(len(cross_section_energy)):
            wavelength = self.energy_to_wavelength(cross_section_energy[i], mass)
            cross_section_wavelength.append(wavelength)
        cross_area_spl = self.cross_section_interpolate(cross_section_sigma, cross_section_energy, mass)
        wavelength_new = np.linspace(min(cross_section_wavelength), max(cross_section_wavelength), 300)
        print('wave', wavelength_new)
        cross_area = cross_area_spl(wavelength_new)

        fac_total = []
        for i in range(len(cross_area)):    
            fac = np.exp( - n * cross_area[i] * distance)
            fac_total.append(fac)

        fac_points = []
        for i in range(len(cross_section_wavelength)):
            fac = np.exp( - n * cross_section_sigma[i] * distance)
            fac_points.append(fac)  

        plt.plot(wavelength_new, fac_total, label = 'interp')
        plt.plot(cross_section_wavelength, fac_points, 'o')
        plt.show()
        return fac 
    
Eu = attenuation_fac(1) 
cross_section_energy = np.linspace(1, 10, num=30)
cross_section_sigma = np.cos(-cross_section_energy**2 / 9.0)
Eu.factor(2,10, 1, cross_section_sigma, cross_section_energy)
