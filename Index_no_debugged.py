import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.spatial import cKDTree
import Structural_fac

# the Xray emission curve 

class index: 
    def __init__(self, image, sample_detetor, cell_len_a, cell_len_b, cell_len_c):
        # image is a string of directory: e.f. "/Users/vivian/Desktop/Research Intern 
        # Oxford/Code/Structure factor/EuMnSb2_untwined2.tif"
        self.image = plt.imread(image)
        self.sample_detector = sample_detetor
        self.cell_len_a = cell_len_a
        self.cell_len_b = cell_len_b
        self.cell_len_c = cell_len_c

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
        pos_rel_center = length * [np.sin(normal_vec[1] / length), np.sin(normal_vec[2] / length)]
        pos = [center[0] + pos_rel_center[0], center[1] + pos_rel_center[1]]
        return pos 
    
    # generate theoretical miller index position map 
    def hkl_map(self, axis, sample_detector_distance, center, hkl_range, domain): 
        # domain: the detecting range of the detector 
        # hkl_range: the range of hkl values
        pos_total = []
        index_total = []
        for h in range(hkl_range[0]):
            for k in range(hkl_range[1]):
                for l in range(hkl_range[2]): 
                    pos = self.hkl_single(axis, sample_detector_distance, h, k, l, center)
                    index = [h, k, l]
                    if abs(pos[0]) < domain[0] and abs(pos[1]) < domain[1]: 
                        pos_total.append(pos)
                        index_total.append(index)

        pos_total = np.array(pos_total)
        index_total = np.array(index_total)
        # plt.matshow(pos)
        # plt.show()
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
    

    def Xray(self):
        # generate a function of Xray intensity versus energy from references 
        Xray_func = 0
        return Xray_func

# cannot do this 
# should be an intergration of all the wavelength 
# need to write he math function 
    # def intensity_singel_wavelength(self, wavelength, Xray_intensity_distribution, structural_val):
    #     h = 6.62607015E-34
    #     energy = h * 1/wavelength
    #     Xray_intensity = Xray_intensity_distribution(energy)
    #     intensity = Xray_intensity * abs(structural_val)
    #     return intensity
    
    # def intensity_total(self, Xray_intensity_distribution, range, )

index1 = index()
distance = index1.dis_pix([1,1], [2, 2])
print(distance)