import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
import Structural_fac

# the Xray emission curve 
def Xray(a): 
    if a >= 12 and a <= 50: 
        if a >= 17.35 and a <= 17.45: 
            Xray_func = lambda a : 4.56087543e-03 * np.exp(-(a - 1.74000000e+01) ** 2 / (2 * (5.80508485e-02 ** 2))) 

        if a >= 19.55 and a <= 19.65: 
              Xray_func = lambda a: 4.78665415e-03 * np.exp(-(a - 1.96000000e+01) ** 2 / (2 * ( 5.61550931e-02 ** 2)))
        else: 
            Xray_func = lambda a: -2E-09* a ** 5 + 4E-07* a ** 4 - 2E-05* a ** 3 + 0.0006* a ** 2 - 0.0073* a + 0.0299
    return Xray_func(a)

print(Xray(19.6))