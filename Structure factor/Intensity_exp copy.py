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
    return intensity_total


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

# This part is to plot out the bright blogs after manual selection. Commenting these out 
# means plotting all the bright blogs detected through the algorithm. 
'''
intensity_list = [[274.        , 174.        ,   7.07106781],    
       [273.        , 697.        ,   7.07106781],
       [170.        , 329.        ,   7.07106781],
       [378.        , 539.        ,   7.07106781],
       [170.        , 539.        ,   7.07106781],
       [377.        , 332.        ,   7.07106781],
       [531.        , 435.        ,   7.07106781],
       [351.        , 434.        ,   7.07106781],
       [ 91.        , 250.        ,   7.07106781],
       [454.        , 617.        ,   7.07106781],
       [454.        , 256.        ,   7.07106781],
       [135.        , 433.        ,   7.07106781],
       [411.        , 434.        ,   7.07106781],
       [477.        , 435.        ,   7.07106781],
       [ 89.        , 618.        ,   7.07106781],
       [273.        , 640.        ,   7.07106781],
       [163.        , 322.        ,   7.07106781],
       [ 68.        , 433.        ,   7.07106781],
       [161.        , 547.        ,   7.07106781],
       [275.        , 230.        ,   7.07106781],
       [274.        , 572.        ,   7.07106781],
       [275.        , 296.        ,   7.07106781],
       [275.        , 215.        ,   7.07106781],
       [455.        , 527.        ,   7.07106781],
       [369.        , 434.        ,   7.07106781],
       [274.        , 551.        ,   7.07106781],
       [272.        , 655.        ,   7.07106781],
       [355.        , 517.        ,   7.07106781],
       [103.        ,  90.        ,   7.07106781],
       [275.        , 359.        ,   7.07106781],
       [196.        , 433.        ,   7.07106781],
       [274.        , 318.        ,   7.07106781],
       [390.        , 434.        ,   7.07106781],
       [356.        , 353.        ,   7.07106781],
       [ 45.        , 318.        ,   7.07106781],
       [159.        , 433.        ,   7.07106781],
       [540.        , 303.        ,   7.07106781],
       [454.        , 345.        ,   7.07106781],
       [275.        , 511.        ,   7.07106781],
       [274.        , 339.        ,   7.07106781],
       [273.        , 626.        ,   7.07106781],
       [275.        , 266.        ,   7.07106781],
       [274.        , 243.        ,   7.07106781],
       [273.        , 604.        ,   7.07106781],
       [274.        , 587.        ,   7.07106781],
       [275.        , 282.        ,   7.07106781],
       [275.        , 113.        ,   7.07106781],
       [180.        , 434.        ,   7.07106781],
       [123.        , 432.        ,   7.07106781],
       [492.        , 435.        ,   7.07106781],
       [ 82.        , 433.        ,   7.07106781],
       [463.        , 435.        ,   7.07106781],
       [426.        , 434.        ,   7.07106781],
       [ 51.        , 433.        ,   7.07106781],
       [443.        , 433.        ,   7.07106781],
       [158.        , 665.        ,   7.07106781],
       [442.        , 777.        ,   7.07106781],
       [ 63.        , 222.        ,   7.07106781],
       [407.        , 706.        ,   7.07106781],
       [160.        , 205.        ,   7.07106781],
       [ 89.        , 526.        ,   7.07106781],
       [386.        , 547.        ,   7.07106781],
       [ 60.        , 646.        ,   7.07106781],
       [481.        , 643.        ,   7.07106781],
       [434.        , 596.        ,   7.07106781],
       [272.        , 760.        ,   7.07106781],
       [386.        , 323.        ,   7.07106781],
       [110.        , 597.        ,   7.07106781],
       [498.        , 662.        ,   7.07106781],
       [181.        , 619.        ,   7.07106781],
       [ 42.        , 664.        ,   7.07106781],
       [ 46.        , 203.        ,   7.07106781],
       [418.        , 580.        ,   7.07106781]]
blobs_log = np.array(intensity_list)
'''

blobs_list = [blobs_log, blobs_dog, blobs_doh]

# remove the edges from blobs_list
blobs_list = remove_edge(blobs_list, image_gray, 10)
# remove the circle from blobs_list
blobs_list = remove_center(blobs_list, image_gray, [430, 280], 50)
# manually remove a few points: 
blobs_list = remove_manual(blobs_list, image_gray, [[819, 163], [88.7, 23.9], [730, 46]])

colors = ['yellow', 'lime', 'red']
titles = ['Laplacian of Gaussian', 'Difference of Gaussian', 'Determinant of Hessian']
sequence = zip(blobs_list, colors, titles )

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
intensity_list = [blobs_log]
I_total = intensity(intensity_list, Iarray)
print(I_total)


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