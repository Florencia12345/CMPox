import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from scipy import optimize 
from mpl_toolkits import mplot3d
# detect blob 
import skimage
from skimage import data 
from skimage.feature import blob_dog, blob_log, blob_doh
# from skimage.color import rgb2gray

# detect circle edges 
from skimage import data, color
from skimage.transform import hough_circle, hough_circle_peaks
from skimage.feature import canny
from skimage.draw import circle_perimeter
from skimage.util import img_as_ubyte
from skimage.transform import hough_ellipse
from skimage.draw import ellipse_perimeter

import matplotlib.pyplot as plt

from skimage import data, color, img_as_ubyte
from skimage.feature import canny
# read params 
# return a, b, c
'''
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
        self.location = [[0, 0, 0], [0.5,   0.5,  0.5]]
        self.locationtotal = [[[0.00000 ,   0.00000 ,   0.00000], [0.00000,    0.50000  ,  0.50000 ], [0.50000  ,  0.00000,    0.50000], [0.50000,    0.50000  ,  0.00000]], [[0.50000,    0.00000 ,   0.00000], [0.5, 0.5, 0.5], [0, 0, 0.5], [0, 0.5, 0]]]
        
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
        self.Eu2 = [self.Eu_a2, self.Eu_b2, self.Eu_c2]

        self.Mn_a = [24.5148, 14.8058, 13.0799, 2.72477]
        self.Mn_b = [2.57255, 0.3493, 15.1280, 146.103]
        self.Mn_c = 7.85731
        self.Mn_a2 = [11.709542 , 1.733414 ,  2.673141,   2.023368,   7.003180]
        self.Mn_b2 = [5.597120,   0.017800,  21.788419 , 89.517915,   0.383054]
        self.Mn_c2 = [-3.750000]
        self.Mn2 = [self.Mn_a2, self.Mn_b2, self.Mn_c2]

        self.Sb_a = [24.5148, 14.8058, 13.0799, 2.72477]
        self.Sb_b = [2.57255, 0.3493, 15.1280, 146.103]
        self.Sb_c = 7.85731
        self.Sb_a2 = [5.394956 , 6.549570,  19.650681,  1.827820 , 17.867833]
        self.Sb_b2 = [33.326523 ,  0.030974,   5.564929,  87.130965,   0.523992]
        self.Sb_c2 = [-0.290506]
        self.Sb2 = [self.Sb_a2, self.Sb_b2, self.Sb_c2]

        self.Cu_a2 = [14.014192 ,  4.784577 ,  5.056806 ,  1.457971,   6.932996]
        self.Cu_b2 = [3.738280 ,  0.003744,  13.034982,  72.554793,   0.265666]
        self.Cu_c2 = [ -3.254477]
        self.Cu2 = [self.Cu_a2, self.Cu_b2, self.Cu_c2]

        self.Pd_a2 = [6.121511,   4.784063 , 16.631683 ,  4.318258,  13.246773]
        self.Pd_b2 = [0.062549 ,  0.784031  , 8.751391 , 34.489983  , 0.784031]
        self.Pd_c2 = [ 0.883099]
        self.Pd2 = [self.Pd_a2, self.Pd_b2, self.Pd_c2]

        self.Cs_a2 = [17.418675,   8.314444 , 10.323193  , 1.383834 , 19.876252]
        self.Cs_b2 = [0.399828,   0.016872,  25.605828, 233.339674 ,  3.826915]
        self.Cs_c2 = [-2.322802]
        self.Cs2 = [self.Cs_a2, self.Cs_b2, self.Cs_c2]

        self.Cl_a2 = [1.446071,   6.870609,  6.151801,   1.750347,   0.634168]   
        self.Cl_b2 = [0.052357,   1.193165,  18.343416 , 46.398394  , 0.401005]
        self.Cl_c2 = [0.146773]
        self.Cl2 = [self.Cl_a2, self.Cl_b2, self.Cl_c2]


# calculate the distance between planes 
    def distance(self, h, k, l): 
        a_recip = 2*np.pi/ self.cell_len_a
        b_recip = 2*np.pi/ self.cell_len_b
        c_recip = 2*np.pi/ self.cell_len_c
        distance = 2*np.pi /(np.sqrt((h * a_recip)* (h * a_recip) + (k * b_recip) * (k * b_recip) + (l * c_recip)* (l * c_recip)))
        return distance

# calculate the brag diffraction angle
    def bragangle(self, wavelength, h, k, l):
        d = self.distance(h, k, l)
        print("distance", d)
        theta = np.arcsin(wavelength/(2*d))
        angle = 2 * theta / (2 * np.pi) * 360
        return angle 
    
    def scat_fac_f0(self, wavelength, element, h, k, l):
        get_ele = eval("self." + element)
        ele_a = get_ele[0]
        ele_b = get_ele[1]
        ele_c = get_ele[2]
        iter = len(get_ele[1])
        theta = self.bragangle(wavelength,h, k, l)
        s = np.sin(theta)/wavelength
        i = 0
        scat_fac = 0
        while i < iter:
            sum = ele_a[i] * np.exp(0 - ele_b[i] * s * s) 
            scat_fac = scat_fac + sum
            # print("in the loop", scat_fac)
            i +=1 
        scat_fac = scat_fac + ele_c[0]
        return scat_fac
    

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

        f = [-11.3+ 7.55621E+00j + f[0],  -11.3+  6.67527E-01j + f[1]]
        print(f)
        for m in range(len(self.locationtotal)): 
            print(m)
            for i in range(len(self.locationtotal[m])): 
                print(i)  
                sum = f[m] * np.exp( 1j * 2 * np.pi * (h * self.locationtotal[m][i][0] + 
                                                       k * self.locationtotal[m][i][1] + 
                                                       l * self.locationtotal[m][i][2]))
                print("location", self.locationtotal[m][i])
                print('sum', sum)
                print('f', f[m])
                S = S + sum 
                print(S)
        return S 


c= crystal("pm", 7.06201,  7.06201, 7.06201, 90, 90, 90)
print(c.distance(0, 0, 2))
# print(c.pnma()[3])
print(c.structure_fac(1.5, ['Cs2', 'Cl2'], 2, 8, 2))

# b = crystal("pnma", 2.98800, 2.98800, 2.98800, 90, 90, 90)
# # print(b.distance(1, 0, 0))
# print(b.structure_fac(1.5, ['Cu2', 'Pd2'], 1, 1, 1))

'''

# input TIFF files
I = plt.imread('/Users/vivian/Desktop/Research Intern Oxford/Code/Structure factor/EuMnSb2_untwined2.tif')
Iarray = np.array(I)
print(Iarray.shape)
Iarray_log = np.log(Iarray)

Iarray_inverse = (Iarray)
Iarray_inverse_crop = Iarray_inverse[240:350, 400:500]
image_rgb  = Iarray_inverse_crop
image_gray = Iarray_inverse_crop
# Iarray_inverse_crop = Iarray_inverse_crop/Iarray_inverse_crop.max()
ave_circ = sum(sum(Iarray_inverse_crop))/(Iarray_inverse_crop.shape[0] * Iarray_inverse_crop.shape[1])

# Iarray_inverse_crop[np.where((Iarray_inverse_crop > ave_circ))] = 1
Iarray_inverse_crop[np.where((Iarray_inverse_crop <= ave_circ))] = 0.2

# center of the detector edeg 
edges = canny(Iarray, sigma=1, low_threshold=0, high_threshold=133)

result = hough_ellipse(edges, accuracy=20, threshold=250, min_size=100, max_size=120)
result.sort(order='accumulator')

# Estimated parameters for the ellipse
best = list(result[-1])
yc, xc, a, b = (int(round(x)) for x in best[1:5])
orientation = best[5]

# Draw the ellipse on the original image
cy, cx = ellipse_perimeter(yc, xc, a, b, orientation)
image_rgb[cy, cx] = (0, 0, 255)
# Draw the edge (white) and the resulting ellipse (red)
edges = color.gray2rgb(img_as_ubyte(edges))
edges[cy, cx] = (250, 0, 0)

fig2, (ax1, ax2) = plt.subplots(
    ncols=2, nrows=1, figsize=(8, 4), sharex=True, sharey=True
)

ax1.set_title('Original picture')
ax1.imshow(image_rgb)

ax2.set_title('Edge (white) and result (red)')
ax2.imshow(edges)

plt.show()

'''
# Detect two radii
hough_radii = np.arange(20, 200, 1)
hough_res = hough_circle(edges, hough_radii)

# Select the most prominent 1 circles
accums, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=3)

# Draw them
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 4))
image = color.gray2rgb(Iarray_inverse_crop)
for center_y, center_x, radius in zip(cy, cx, radii):
    circy, circx = circle_perimeter(center_y, center_x, radius, shape=image.shape)
    image[circy, circx] = (220, 20, 20)

ax.imshow(Iarray_inverse_crop, cmap=plt.cm.gray)
plt.show()
'''