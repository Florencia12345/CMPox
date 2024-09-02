import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.spatial import cKDTree
import Structural_fac
import Xray_emit_curve

# This class is responsible for calculating the location of bright peaks based on their indices, e.g., [1, 0, 14].
# Note: The function to calculate the theoretical intensity on Laue patterns is still under development.
class index: 
    def __init__(self, sample_detetor, cell_len_a, cell_len_b, cell_len_c):
        # image is a string of directory: e.f. "/Users/vivian/Desktop/Research Intern 
        # Oxford/Code/Structure factor/EuMnSb2_untwined2.tif"
        # self.image = plt.imread(image)      if used, then need to add another variable as "image"
        self.sample_detector = sample_detetor
        self.cell_len_a = cell_len_a
        self.cell_len_b = cell_len_b
        self.cell_len_c = cell_len_c
        self.wavelength_range = [0, 1.24]
        self.Xray = Xray_emit_curve.Xray

    def center_detector(self, center): 
    # find the center of the diffraction peak from the diffraction pattern 
        center1 = center
        return center1

    def angle_pix(self,center, pos):
        center = self.center_detector(center)
        center = np.array(center)
        pos = np.array(pos)
        distance1 = np.sqrt(np.dot((center - pos), (center-pos)))
        distance2 = self.sample_detector
        alpha = np.atan(distance1/distance2)
        theta = (np.pi - alpha)/2 
        return theta 
    
    '''
    def hkl_map2(self, wavelength, h, k, l, center, pos):
        # axis: the beam incoming axis. e.f. EuMnSb2: incoming beam is the a axis 
        theta = self.angle_pix(center, pos)
        d = wavelength / (2 * np.sin(theta))
        return 
    '''
    # @staticmethod
    def Xray(self, a): 
        if a >= 12 and a <= 50: 
            if a >= 17.35 and a <= 17.45: 
                Xray_func = lambda a : 4.56087543e-03 * np.exp(-(a - 1.74000000e+01) ** 2 / (2 * (5.80508485e-02 ** 2))) 

            if a >= 19.55 and a <= 19.65: 
                Xray_func = lambda a: 4.78665415e-03 * np.exp(-(a - 1.96000000e+01) ** 2 / (2 * ( 5.61550931e-02 ** 2)))
            else: 
                Xray_func = lambda a: -2E-09* a ** 5 + 4E-07* a ** 4 - 2E-05* a ** 3 + 0.0006* a ** 2 - 0.0073* a + 0.0299
            return Xray_func(a)
        else: 
            print('out of the range')
            return 0 
        
    # calculate a single hkl miller index position 
    def hkl_single(self, axis, sample_detector_distance, h, k, l, center):
        # center: the center of diffraaction pattern
        # axis: the imcoming beam axis in a, b, or c 

        a_recip = 2*np.pi/ self.cell_len_a
        b_recip = 2*np.pi/ self.cell_len_b
        c_recip = 2*np.pi/ self.cell_len_c
        # distance = 2*np.pi /(np.sqrt((h * a_recip)**2 + (k * b_recip)**2 + (l * c_recip)**2))
        # needs to be a 3D projection of angle to a 2D plane 
        normal_vec = [h * a_recip, k * b_recip, l * c_recip]
        if axis == "a":
            beam_vec = np.array([1, 0, 0])
            direction = [k, l]
        if axis == "b":
            beam_vec = np.array([0, 1, 0])
            direction = [h, l]
        if axis == "c":
            beam_vec = np.array([0, 0, 1])
            direction = [h, k]
        # calculate the angle between beam direction and the nomral vector of the plane
        alpha = np.arccos(np.dot(beam_vec, normal_vec)/(np.linalg.norm(beam_vec) * np.linalg.norm(normal_vec)))
        angle = 2 * alpha
        length = sample_detector_distance * np.tan(angle)
        pos_rel_center = [length * np.sin(np.arctan(normal_vec[1] / normal_vec[2])), length * np.sin(np.arctan(normal_vec[2] / normal_vec[1]))]
        pos2pixel = [pos_rel_center[0] * 255/46, pos_rel_center[1] * 255/46]
        # pos_rel_center = length * [np.sin(normal_vec[1] / length), np.sin(normal_vec[2] / length)]
        pos = [center[0] + pos2pixel[0], center[1] + pos2pixel[1]]
        return pos, alpha
    
    # generate theoretical miller index position map 
    def hkl_map(self, axis, sample_detector_distance, center, hkl_range, domain): 
        # domain: the detecting range of the detector 
        # hkl_range: the range of hkl values
        pos_total = []
        index_total = []
        for h in range(hkl_range[0]):
            for k in range(hkl_range[1]):
                for l in range(hkl_range[2]): 
                    pos, alpha = self.hkl_single(axis, sample_detector_distance, h, k, l, center)
                    index = [h, k, l]
                    if abs(pos[0]) < domain[0] and abs(pos[1]) < domain[1]: 
                        pos_total.append(pos)
                        index_total.append(index)
        pos_total = np.array(pos_total)
        index_total = np.array(index_total)
        return pos_total, index_total
    
 
    # compare the experimental data with the theoretical data 
    # the exp_data needs to be pre-treated so that it captures the same bright peaks as the theo_data. 
    # i.e. length should be the same, and points should be the same through manual supervision. 
    def match(exp_data, theo_data, index_total):
        # theo_data is the theoretically calculated positions of the miller index from function hkl_map
        # cKD maps 
        tree = cKDTree(exp_data)
        distances, indices = tree.query(theo_data)
        # create a mapping dictionary where the keys are the indixes of the theoretical points 
        # and the values are the indices of the nearest experimental points
        mapping = {i:indices[i] for i in range(len(theo_data))}
        miller_mapping = {index_total[i]: exp_data[indices[i]] for i in range(len(theo_data))}
        return mapping, miller_mapping
    
    # needs further amendation
    def intensity_hkl(self, element, h, k, l, axis, sample_detector_distance, center):
        # wavelength = self.wavelength_range
        crystal_hkl = Structural_fac.crystal("pnma", 22.567, 4.371, 4.409, 90, 90, 90)
        distance = crystal_hkl.distance(h, k, l)
        # compute the brag angle
        pos, alpha = self.hkl_single(axis, sample_detector_distance, h, k, l, center)
        flag = 1
        n = 1
        intensity = 0
        print('alpha', alpha)
        while flag: 
            wavelength = np.sin(alpha) * 2 * distance * n 
            print('wavelength A', wavelength)
            energy = 6.62607015E-34 * 3E8 / (wavelength * 1E-10)
            keV = energy / 1.602E-16
            print('energy', keV)

            if wavelength >= min(self.wavelength_range) and wavelength <= max(self.wavelength_range):
                structure_fac = crystal_hkl.structure_fac(wavelength, element, h, k, l)
                print(structure_fac)
                ratio = self.Xray(keV)
                intensity = intensity + abs(structure_fac) ** 2 * ratio
                print('intensity', intensity)
            if wavelength >= max(self.wavelength_range): 
                flag = 0
            n += 1 
        return intensity 
    
    def intensity_list(self, element, hkllist, axis, sample_detector_distance, center):
        intlist = []
        for i in range(len(hkllist)):
            print(i)
            h = hkllist[i][0]
            k = hkllist[i][1]
            l = hkllist[i][2]
            Iy = self.intensity_hkl(element, h, k, l, axis, sample_detector_distance, center)
            intlist.append(Iy)
            print(Iy)
        return intlist
    
index1 = index(1, 22.567, 4.371, 4.409)
a, alpha = index1.hkl_single('a', 55, 14, -1, .0000000001, [432, 274])
print('a', a)
    
intlist = [[17, -1, 1],
[16, 0.000000000000000001, -1],
[18, -1, 1],
[17, 0.000000000000000001, -1],
[18, 0.000000000000000001, -1],
[20, -1, 1],
[22, 0.000000000000000001, -1],
[24, 0.000000000000000001, -1],
[28, 0.000000000000000001, -1],
[30, -1, 1],
[32, -1, 1],
[32, 0.000000000000000001, -1],
[16, 1, 0.000000000000000001],
[12, 1, 0.000000000000000001],
[20, 1, 0.000000000000000001],
[18, 1, 0.000000000000000001],
[17, 1, 0.000000000000000001],
[14, 1, 0.000000000000000001],
[14, -1, 0.000000000000000001],
[18, -1, 0.000000000000000001],
[28, -1, 0.000000000000000001],
[32, -1, 0.000000000000000001],
[28, 1,0.000000000000000001],
[24, 1, 0.000000000000000001],
[22, 1, 0.000000000000000001],
[12, -1, 0.000000000000000001],
[16, -1, 0.000000000000000001],
[17, -1, 0.000000000000000001],
[20, -1, 0.000000000000000001],
[22, -1, 0.000000000000000001],
[24, -1, 0.000000000000000001],
[38, -1, 0.000000000000000001],
[32, 1, 0.000000000000000001],
[38, 0.000000000000000001, 1],
[32, 0.000000000000000001, 1],
[28, 0.000000000000000001, 1],
[24, 0.000000000000000001, 1],
[24, 1, 1],
[22, 0.000000000000000001, 1],
[20, 0.000000000000000001, 1],
[20, 1, 1],
[18, 0.000000000000000001, 1],
[17, 0.000000000000000001, 1],
[18, 1, 1],
[16, 0.000000000000000001, 1],
[17, 1, 1],
[14, 0.000000000000000001, 1]]

# intensity = index1.intensity_hkl(['Eu2', 'Mn2', 'Sb2', 'Sb2'], 24, 0.0000000001, 1, 'a', 55, [436, 286])
# inttotal = index1.intensity_list(['Eu2', 'Mn2', 'Sb2', 'Sb2'], intlist, 'a', 55, [436, 286])
# distance = index1.dis_pix([1,1], [2, 2])
