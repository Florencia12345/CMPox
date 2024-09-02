import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

# df = pd.read_excel('/Users/vivian/Desktop/Research Intern Oxford/Squid data/EuAuAs/EuAuAs_Hc_Tscan.xlsx') 
df = pd.read_excel('/Users/vivian/Desktop/Research Intern Oxford/Squid data/EuAuAs/EuAuAs_Hc_zfc_Tscan.xlsx') 
# df = pd.read_excel('/Users/vivian/Desktop/Research Intern Oxford/Squid data/EuMnSb2/EuMnSb2_Tscan_zfc_Ha.xlsx') 

temp = df['Temperature (K)']
moment = df['Moment (emu)']
field_h = df['Magnetic Field (Oe)']

temp = np.array(temp)
moment = np.array(moment)
field_h = np.array(field_h)

# normalization, convert suscept_SI into suscept(emu/mol)
# na = 6.022E23
# mu_0 = 1.25663706E-6
mass = 62.2E-3
# mass = 2.7E-3
mmol = 151.964 + 196.96657 + 74.9216
# mmol = 151.964 + 54.94 + 121.7 * 2
n = mass/mmol 
moment = moment / n
suscept = moment / field_h

print(suscept)

def Moment_vs_Temp(temp, moment):
    plt.plot(temp, moment)
    plt.title('Moment_vs_Temp')
    plt.xlabel("Temperature(K)")
    plt.ylabel("Induced Moment(emu)")
    plt.show()


def Sus_inverse_vs_Temp(suscept, temp):
    sus_inverse = [1/ x for x in suscept]
    plt.plot(temp, sus_inverse)
    plt.scatter(temp, sus_inverse)
    plt.title('1/X_vs_Temp')
    plt.xlabel("Temperature(K)")
    plt.ylabel("1/X")
    plt.show()

def cutoff(suscept, temp): 
    sus_inverse = [1/ x for x in suscept]
    length = len(sus_inverse)
    linear_sus = sus_inverse[(round(length/3) - 1):]
    linear_temp = temp[(round(length/3) - 1):]
    coef = np.polyfit(linear_temp, linear_sus, 1)
    print(coef)
    poly1d_fn = np.poly1d(coef) 
    plt.plot(temp, sus_inverse, 'yo', linear_temp, poly1d_fn(linear_temp), '--k')
    plt.title('1/X_vs_Temp')
    text = '1/X = ' + str(round(coef[0], 10)) + '*' + 'T' + str(round(coef[1], 6))
    plt.text(max(temp)/2, max(sus_inverse)/3*2, text)
    plt.xlabel("Temperature(K)")
    plt.ylabel("1/X")

    plt.show()
    return coef

coef = cutoff(suscept, temp)
C_prime = 1/ coef[0]
T_c = -coef[1] * C_prime

print(C_prime, T_c)

