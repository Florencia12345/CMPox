import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

class curie_weiss: 
    def __init__(self, filename, mass, mmol) -> None:
        self.df = pd.read_excel(filename)
        self.temp = np.array(self.df['Temperature (K)'])
        self.moment = np.array(self.df['Moment (emu)'])
        self.field_h = np.array(self.df['Magnetic Field (Oe)'])
        self.mass = mass
        self.mmol = mmol
        self.n = mass/mmol
        self.suscept = self.moment/ self.n / self.field_h

    def Moment_vs_Temp(self):
        plt.plot(self.temp, self.moment)
        plt.title('Moment_vs_Temp')
        plt.xlabel("Temperature(K)")
        plt.ylabel("Induced Moment(emu)")
        plt.show()

    def Sus_inverse_vs_Temp(self):
        sus_inverse = [1/ x for x in self.suscept]
        plt.plot(self.temp, sus_inverse)
        plt.scatter(self.temp, sus_inverse)
        plt.title('1/X_vs_Temp')
        plt.xlabel("Temperature(K)")
        plt.ylabel("1/X")
        plt.show()

    def cutoff(self): 
        sus_inverse = [1/ x for x in self.suscept]
        length = len(sus_inverse)
        linear_sus = sus_inverse[(round(length/3) - 1):]
        linear_temp = self.temp[(round(length/3) - 1):]
        coef = np.polyfit(linear_temp, linear_sus, 1)
        print(coef)
        poly1d_fn = np.poly1d(coef) 
        plt.plot(self.temp, sus_inverse, 'yo', linear_temp, poly1d_fn(linear_temp), '--k')
        plt.title('1/X_vs_Temp')
        text = '1/X = ' + str(round(coef[0], 10)) + '*' + 'T' + str(round(coef[1], 6))
        plt.text(max(self.temp)/2, max(sus_inverse)/3*2, text)
        plt.xlabel("Temperature(K)")
        plt.ylabel("1/X")
        plt.show()
        return coef
    
    def ct(self):
        coef = self.cutoff()
        C_prime = 1/ coef[0]
        T_c = -coef[1] * C_prime
        print(C_prime, T_c)
        return [C_prime, T_c]

curie = curie_weiss('/Users/vivian/Desktop/Research Intern Oxford/Squid data/EuAuAs/EuAuAs_Hc_Tscan.xlsx', 62.2E-3, 151.964 + 196.96657 + 74.9216)
ct = curie.ct()