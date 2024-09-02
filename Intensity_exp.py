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


# Gaussian fit
def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-x)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-y)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

# create mask for the edges 
def edge_mask(image,width, height): 
    size = image.shape 
    frame = np.zeros(size)
    frame[0 : width, :] = image[0 : width, :]
    frame[:,0: height] = image[:,0: height]
    frame[(size[0]-width) : size[0], :] = image[(size[0]-width) : size[0], :]
    frame[:, (size[1] - height) : size[1]] = image[:, (size[1] - height) : size[1]]
    return frame

def remove_edge(blobs_list, image, width):
    shape = image.shape
    lenth = [shape[0], shape[1]]
    pos = [[] for x in range(len(blobs_list))]
    for i in range(len(blobs_list)):
        for ii in range(len(blobs_list[i])):
            if blobs_list[i][ii][1] < width or blobs_list[i][ii][1] > (lenth[1] - width) \
            or blobs_list[i][ii][0] < width or blobs_list[i][ii][0] > (lenth[0] - width):
                pos[i].append(ii)
    for i in range(len(blobs_list)):                
        blobs_list[i] = np.delete(blobs_list[i], pos[i], axis = 0)

    return blobs_list

def remove_center(blobs_list, image, center, radi):
    # fill in the center of the detector and radius, remove the spotted bright blogs
    shape = image.shape
    lenth = [shape[0], shape[1]]
    pos = [[] for x in range(len(blobs_list))]
    for i in range(len(blobs_list)):
        for ii in range(len(blobs_list[i])):
            if ((blobs_list[i][ii][1]- center[0])**2 + (blobs_list[i][ii][0]-center[1])**2) <  radi **2:
                pos[i].append(ii)
                print(pos)

    for i in range(len(blobs_list)):                
        blobs_list[i] = np.delete(blobs_list[i], pos[i], axis = 0)
    
    return blobs_list

def remove_manual(blobs_list, image, pos):
    for iii in range(len(pos)): 
        pos = pos[iii]
        position = [[] for x in range(len(blobs_list))]
        for i in range(len(blobs_list)):
            for ii in range(len(blobs_list[i])):
                if ((blobs_list[i][ii][1]- pos[0])**2 + (blobs_list[i][ii][0]-pos[1])**2) <  100:
                    print(blobs_list[i][ii])
                    print(i)
                    print(ii)
                    position[i].append(ii)

        print(position)

        for i in range(len(blobs_list)):                
            blobs_list[i] = np.delete(blobs_list[i], position[i], axis = 0)
        return blobs_list

# calculate the intensity 
def intensity(blobs_list,image): 
    intensity_total = []
    for i in range(len(blobs_list)):
    # get the first array
        position = blobs_list[i] #i: 0,1,2 for three different detection,
        for j in range(len(blobs_list[i])):
            single_pos = position[j] #x, y and radius 
            if abs(single_pos[0]-image.shape[0])>single_pos[2] and single_pos[0] > single_pos[2] \
            and abs(single_pos[1]-image.shape[1])>single_pos[2] and single_pos[1] > single_pos[2]:
                domain = [[int(single_pos[0]- round(single_pos[2])), int(single_pos[0] + round(single_pos[2]))], 
                        [int(single_pos[1] - round(single_pos[2])), int(single_pos[1] + round(single_pos[2]))]]
                image_crop = image[domain[0][0]:domain[0][1], domain[1][0]:domain[1][1]]

            # do a 2d gaussian fit --- gaussian fit doesn't look good
            # params = fitgaussian(image_crop)
            # fit = gaussian(*params)
            # fit_value = fit(*np.indices(image_crop.shape))
            # plt.matshow(fit_value)
            # plt.show()
            
                # calculate the average of edge(which is similar to the background) and subtract it 
                edge = np.zeros(image_crop.shape)
                edge[0, :(edge.shape[0]-1)] = image_crop[0, :(edge.shape[0]-1)]
                edge[(edge.shape[0]-1), :(edge.shape[0]-1)] = image_crop[(edge.shape[0]-1), :(edge.shape[0]-1)]
                edge[:(edge.shape[0]-1), 0] = image_crop[: (edge.shape[0]-1), 0]
                edge[:(edge.shape[0]), (edge.shape[0]-1)] = image_crop[:(edge.shape[0]), (edge.shape[0]-1)]
                ave = sum(sum(edge) )/ ((image_crop.shape[0]+image_crop.shape[1]-2)*2)
                # plt.matshow(edge)
                # plt.show()
                image_remove_err = image_crop - np.ones(image_crop.shape) * ave 
                # plt.matshow(image_remove_err)
                # plt.show()
                intensity = sum(sum(image_remove_err))
                coord = [i, j, intensity]
                intensity_total.append(coord)
                # print(intensity_total)
            else:
                print("Too close to the boundary")
    return image_crop, blobs_list



# input TIFF files
I = plt.imread('/Users/vivian/Desktop/Research Intern Oxford/Code/Structure factor/EuMnSb2_untwined2.tif')
Iarray = np.array(I)
Iarray_log = np.log(Iarray)

# Input TIFF files 
I_ski = skimage.io.imread('/Users/vivian/Desktop/Research Intern Oxford/Code/Structure factor/EuMnSb2_untwined2.tif')
image_gray = I_ski

# remove the Gaussian background distribution 
params = fitgaussian(Iarray)
fit = gaussian(*params)
image_gray = I_ski - fit(*np.indices(Iarray.shape))

# it doesnt work 
# edge_mask1 = edge_mask(image_gray, 2, 2)

print("image_gray", image_gray)
image_gray_crop = (image_gray*1000)

# laplacian of Gaussian and # Compute radii in the 3rd column.
blobs_log = blob_log(image_gray, min_sigma = 5, max_sigma=10, num_sigma=5, threshold=1)
blobs_log[:, 2] = blobs_log[:, 2] * np.sqrt(2)
# Difference of Gaussian
blobs_dog = blob_dog(image_gray, min_sigma = 3.3, max_sigma=10, threshold=1)
blobs_dog[:, 2] = blobs_dog[:, 2] * np.sqrt(2)
# determinant of hessian
blobs_doh = blob_doh(image_gray, min_sigma = 5, max_sigma=10, threshold=1)
blobs_list = [blobs_log, blobs_dog, blobs_doh]

# remove the edges from blobs_list
blobs_list = remove_edge(blobs_list, image_gray, 10)
# remove the circle from blobs_list
blobs_list = remove_center(blobs_list, image_gray, [430, 280], 50)
# manually remove a few points: 
blobs_list = remove_manual(blobs_list, image_gray, [[819, 163], [88.7, 23.9], [730, 46]])



colors = ['yellow', 'lime', 'red']
titles = ['Laplacian of Gaussian', 'Difference of Gaussian', 'Determinant of Hessian']
sequence = zip(blobs_list, colors, titles)

# plot 
fig, axes = plt.subplots(1, 3, figsize=(9, 3), sharex=True, sharey=True)
ax = axes.ravel()
print("blobs_list", blobs_list)
for idx, (blobs, color, title) in enumerate(sequence):
    ax[idx].set_title(title)
    ax[idx].imshow(image_gray)
    for blob in blobs:
        y, x, r = blob
        c = plt.Circle((x, y), r, color=color, linewidth=0.5, fill=False)
        ax[idx].add_patch(c)
    ax[idx].set_axis_off()
plt.tight_layout()
plt.show()

# compute the intensity of all 
[image_crop, intensity] = intensity(blobs_list, Iarray)
















'''
# plot countour 
plt.contour(fit(*np.indices(data.shape)), cmap=plt.cm.copper)
ax = plt.gca()
(height, x, y, width_x, width_y) = params

# plot dat after removing noise 
# data = data -fit(*np.indices(data.shape))
plt.matshow(data)
print(data)

# plot the parameters
plt.text(0.95, 0.05, """
x : %.1f
y : %.1f
width_x : %.1f
width_y : %.1f""" %(x, y, width_x, width_y),
        fontsize=16, horizontalalignment='right',
        verticalalignment='bottom', transform=ax.transAxes)
plt.show()

# # 3D plotting
# y = np.linspace(0, 562, 563)
# x = np.linspace(0, 862, 863)
# X, Y = np.meshgrid(x, y)
# Z = data

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.contour3D(X, Y, Z)
# plt.show()

'''