import numpy as np
# read params 
# return a, b, c

class crystal: 
    def __init__(self, sym, cell_len_a, cell_len_b, cell_len_c, cell_angle_alpha, cell_angle_beta, cell_angle_gamma):
        # define the symmetry and parameters of the crystal
        # e.g. sym = pnma, cell_len_a = 22.567
        self.sym = sym
        self.cell_len_a = cell_len_a
        self.cell_len_b = cell_len_b
        self.cell_len_c = cell_len_c
        self.alpha = cell_angle_alpha
        self.beta = cell_angle_beta
        self.gamma = cell_angle_gamma
        # self.location = [[0, 0, 0], [0.5,   0.5,  0.5]]
        # self.locationtotal = [[[0, 0, 0], [0, 1, 0], [0, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 0], [1, 0, 1], [1, 1, 1]], [[0.5, 0.5, 0.5]]]
        # self.locationtotal = [[[0.5,   0.5,  0.5], [0.5,   -0.5,  0.5], [0.5,   0.5,  -0.5], [-0.5,   0.5,  0.5], [-0.5, -0.5,  -0.5], [0.5, -0.5,  -0.5], [-0.5, 0.5,  -0.5], [-0.5, -0.5,  0.5]], [[0, 0, 0]]]                 
        # self.locationtotal = [[[0, 0, 0]], [[0.5, 0.5, 0.5]]]
                      
        self.location = [[0.38637, 0.25000 , 0.76977 ], [0.24970  ,  0.25000   , 0.27005], 
                         [0.17468  ,  0.25000 ,   0.77005], [0.49905 ,   0.25000 ,   0.27819]]
        
        self.locationtotal = [[[0.38637, 0.25, 0.76977], [0.11363000000000001, 0.75, 1.2697699999999998], [0.88637, 0.25, -0.26976999999999995], [-0.38637, -0.25, -0.76977]],
                              [[0.2497, 0.25, 0.27005], [0.25029999999999997, 0.75, 0.77005], [0.7497, 0.25, 0.22995], [-0.2497, 0.75, -0.27005]],
                                [[0.17468, 0.25, 0.77005], [0.67468, 0.25, -0.27005], [-0.17468, -0.25, -0.77005], [0.32532, 0.75, 1.27005]],
                                [[0.49905, 0.25, 0.27819], [0.0009500000000000064, 0.75, 0.7781899999999999], [0.99905, 0.25, 0.22181], [-0.49905, -0.25, -0.27819]]
        ]

        self.Eu_a = [24.5148, 14.8058, 13.0799, 2.72477]
        self.Eu_b = [2.57255, 0.3493, 15.1280, 146.103]
        self.Eu_c = [7.85731]
        self.Eu1 = [self.Eu_a, self.Eu_b, self.Eu_c]

        self.Eu_a2 = [17.186195, 37.156838,  13.103387 ,  2.707246 , 24.419271]
        self.Eu_b2 = [ 0.261678,  0.001995 , 14.787360, 134.816293 ,  2.581883]
        self.Eu_c2 = [-31.586687]
        self.Eu_bc = [7.220000]
        self.Eu2 = [self.Eu_a2, self.Eu_b2, self.Eu_c2, self.Eu_bc]
        
        self.Mn_a = [24.5148, 14.8058, 13.0799, 2.72477]
        self.Mn_b = [2.57255, 0.3493, 15.1280, 146.103]
        self.Mn_c = 7.85731
        self.Mn_a2 = [11.709542 , 1.733414 ,  2.673141,   2.023368,   7.003180]
        self.Mn_b2 = [5.597120,   0.017800,  21.788419 , 89.517915,   0.383054]
        self.Mn_c2 = [-0.147293]
        self.Mn_bc = [-3.750000]
        self.Mn2 = [self.Mn_a2, self.Mn_b2, self.Mn_c2, self.Mn_bc]

        self.Sb_a = [24.5148, 14.8058, 13.0799, 2.72477]
        self.Sb_b = [2.57255, 0.3493, 15.1280, 146.103]
        self.Sb_c = 7.85731
        self.Sb_a2 = [5.394956 , 6.549570,  19.650681,  1.827820 , 17.867833]
        self.Sb_b2 = [33.326523 ,  0.030974,   5.564929,  87.130965,   0.523992]
        self.Sb_c2 = [-0.290506]
        self.Sb_bc = [5.570000]
        self.Sb2 = [self.Sb_a2, self.Sb_b2, self.Sb_c2, self.Sb_bc]

        self.Cu_a2 = [14.014192 ,  4.784577 ,  5.056806 ,  1.457971,   6.932996]
        self.Cu_b2 = [3.738280 ,  0.003744,  13.034982,  72.554793,   0.265666]
        self.Cu_c2 = [ -3.254477]
        self.Cu2 = [self.Cu_a2, self.Cu_b2, self.Cu_c2]

        self.Pd_a2 = [6.121511,   4.784063 , 16.631683 ,  4.318258,  13.246773]
        self.Pd_b2 = [0.062549 ,  0.784031  , 8.751391 , 34.489983  , 0.784031]
        self.Pd_c2 = [ 0.883099]
        self.Pd2 = [self.Pd_a2, self.Pd_b2, self.Pd_c2]
        


# calculate the distance between planes 
    def distance(self, h, k, l): 
        a_recip = 2*np.pi/ self.cell_len_a
        b_recip = 2*np.pi/ self.cell_len_b
        c_recip = 2*np.pi/ self.cell_len_c
        distance = 2*np.pi /(np.sqrt((h * a_recip)**2 + (k * b_recip)**2 + (l * c_recip)**2))
        return distance

# calculate the brag diffraction angle
    def bragangle(self, wavelength, h, k, l):
        d = self.distance(h, k, l)
        theta = np.arcsin(wavelength*0.00000001/(2*d*0.00000001))
        angle = 2 * theta / (2 * np.pi) * 360
        return angle 
    
    def scat_fac_f0(self, wavelength, element, h, k, l):
        get_ele = eval("self." + element)
        ele_a = get_ele[0]
        ele_b = get_ele[1]
        ele_c = get_ele[2]
        iter = len(get_ele[1])
        theta = self.bragangle(wavelength,h, k, l)/2
        s = np.sin(theta/360 * 2 * np.pi)/(wavelength)
        print("s", s)
        i = 0
        scat_fac = 0
        while i < iter:
            sum = ele_a[i] * np.exp( - ele_b[i] * s**2)
            scat_fac = scat_fac + sum
            i += 1
        scat_fac = scat_fac + ele_c[0]
        return scat_fac
    
    def scat_fac_f1(self, wavelength):
        scat_fac_f1 =0 
        return scat_fac_f1
    
    def scat_fac_f2(self, wavelength):
        scat_fac_f2 = 0
        return scat_fac_f2
    
    def structure_fac(self, wavelength, element, h, k, l):
        # element is a list containing all elements of this material 
        # element = ['Eu2', 'Mn', 'Sb'  'Sb']
        iter = len(self.locationtotal[1])
        f = []
        S = 0
        
        for i in range(len(self.locationtotal)):
            ele = element[i]
            f.append(self.scat_fac_f0(wavelength, ele, h, k, l))
        print(f)

        f = [-6.051638 +12.4478j + f[0] -1.43280E-02, 
             -0.4843259+2.6828j + f[1] -6.24090E-03, 
             -1.3316869999999998+5.62234j + f[2] -1.17200E-02, 
             -1.3316869999999998+5.62234j + f[3] -1.17200E-02]
        # f = [f[0] + f1[0] + f2[0]].....
        # f = [2.27847E+00 + 5.57718E-01j + f[0], 4.62570E-02 + 3.74793E+00j + f[1]]
        # print(f)
        for m in range(len(self.locationtotal)): 
            # print(m)
            for i in range(len(self.locationtotal[m])): 
                # print(i)  
                sum = f[m] * np.exp( 1j * 2 * np.pi * (h * self.locationtotal[m][i][0] + 
                                                       k * self.locationtotal[m][i][1] + 
                                                       l * self.locationtotal[m][i][2]))
                print("location", self.locationtotal[m][i])
                # print('sum', sum)
                # print('f', f[m])
                S = S + sum 
                print(S)
        return S 
    
    def pnma(self): 
        location = []
        for i in range(len(self.location)):
            a = self.location[i]
            group1 = self.location[i]
            group2 = [1/2 - a[0], 1/2 + a[1], 1/2 + a[2]]
            group3 = [a[0], 1/2 - a[1], a[2]]
            group4 = [1/2 + a[0], a[1], 1/2 -a[2]]
            group5 = [-a[0], -a[1], -a[2]]
            group6 = [1/2 + a[0], 1/2 - a[1], 1/2 - a[2]]
            group7 = [-a[0], 1/2 + a[1], -a[2]]
            group8 = [1/2 -a[0], -a[1], 1/2 + a[2]]
            group = [group1, group2, group3, group4, group5, group6, group7, group8]
            location.append(group)

        return location

c= crystal("pnma", 22.567, 4.371, 4.409, 90, 90, 90)
print(c.bragangle(1.5, 4, 0, 0))
# print(c.pnma()[3])
print(c.structure_fac(1.5, ['Eu2', 'Mn2', 'Sb2','Sb2'], 1, 0, 1))

# b = crystal("pnma", 2.98800, 2.98800, 2.98800, 90, 90, 90)
# # print(b.distance(1, 0, 0))
# print(b.structure_fac(1.5, ['Cu2', 'Pd2'], 1, 1, 1))

