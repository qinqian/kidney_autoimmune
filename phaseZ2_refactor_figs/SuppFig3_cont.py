# ###https://jef.works/STalign/notebooks/xenium-heimage-alignment.html

from STalign import STalign
import ants
import scanpy as sc
import numpy as np
from skimage.transform import warp
from skimage.transform import AffineTransform, warp, estimate_transform


# import dependencies
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch
from skimage import io, color, filters, measure, morphology
import numpy as np
from skimage.transform import warp


# make plots bigger
plt.rcParams["figure.figsize"] = (12,10)


def cropped_HE(image_file):

    #V = plt.imread(image_file)

    image = io.imread(image_file)

    # If image has an alpha channel, drop it
    if image.shape[-1] == 4:
        image = image[:, :, :3] 

    # 1. Convert to grayscale
    gray = color.rgb2gray(image) if image.ndim == 3 else image

    # 2. Threshold to separate foreground from background
    thresh = filters.threshold_otsu(gray)
    binary = gray < thresh  # tissue often darker than background

    # Optional: remove small artifacts and fill holes
    binary = morphology.remove_small_objects(binary, min_size=500)
    binary = morphology.remove_small_holes(binary, area_threshold=500)

    # 3. Label connected components
    labeled = measure.label(binary)
    regions = measure.regionprops(labeled)

    # 4. Find the largest component
    largest_region = max(regions, key=lambda r: r.area)

    # 5. Get bounding box (min_row, min_col, max_row, max_col)
    minr, minc, maxr, maxc = largest_region.bbox

    pad_y = int((maxr - minr) * 0.001)
    pad_x = int((maxc - minc) * 0.001)

    h, w = image.shape[:2]

    minr = max(minr - pad_y, 0)
    maxr = min(maxr + pad_y, h)
    minc = max(minc - pad_x, 0)
    maxc = min(maxc + pad_x, w)

    # 6. Crop original image
    cropped = image[minr:maxr, minc:maxc]

    # plot
    fig,ax = plt.subplots(1, 2)
    #ax.imshow(V)
    ax[0].imshow(image)
    ax[1].imshow(cropped)
    Inorm = STalign.normalize(cropped)
    ax[1].imshow(Inorm)    
    fig.savefig("test_HE1.png")
    return Inorm

#image_file = 'BD2025_15_3__2025-07-29_23.04.43.ndpi.png'

image_file = 'BS22_12012A1_small.png'
image = io.imread(image_file)

fig,ax = plt.subplots(5, 1) #, sharey=True)
fig.set_size_inches(8, 32)

cropped = image

ax[0].imshow(image, alpha=1, origin='lower')
#ax[0].invert_yaxis()
#ax[0].set_xlim(1000, 4000)
#ax[0].set_ylim(1000, 2400)
ax[1].imshow(image, alpha=1, origin='lower')
ax[1].set_xlim(2500, 2800)
ax[1].set_ylim(1300, 1600)

def process_xenium_cells(fname):
    # Single cell data to be aligned
    df = pd.read_csv(fname)
    # get cell centroid coordinates
    xM = np.array(df['x_centroid'])
    yM = np.array(df['y_centroid'])

    # rasterize at 30um resolution (assuming positions are in um units) and plot
    XJ,YJ,M,fig = STalign.rasterize(xM, yM, dx=30)
    ax = fig.axes[0]
    ax.invert_yaxis()    
    # plot
    fig,ax = plt.subplots(1, 1)
    ax.scatter(xM,yM,s=1,alpha=0.2)
    return df, XJ, YJ, M, xM, yM

fname = '../data/kidney/20240803__182820__BWH_20240803_skin_Shruti_kidney/output-XETG00150__0018462__BS22_12012A1__20240803__183643/cells.csv.gz'
df, XJ, YJ, M, xM, yM = process_xenium_cells(fname)

I = cropped.transpose(2,0,1)
print(I.shape)
YI = np.array(range(I.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
XI = np.array(range(I.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
extentI = STalign.extent_from_x((YI,XI))

print(M.shape)
J = np.vstack((M, M, M)) # make into 3xNxM
print(J.min())
print(J.max())

# normalize
J = STalign.normalize(J)
print(J.min())
print(J.max())

# double check size of things
print(I.shape)
print(M.shape)
print(J.shape)


# manually make corresponding points
pointsI = np.array([[1000, 200], [1000, 2200], [900, 4500], [50, 6500]])
pointsJ = np.array([[1600, 20], [1500, 2500], [1330, 5300], [200, 8000]])


L,T = STalign.L_T_from_points(pointsI, pointsJ)

pointsJ_trans = L @ pointsI.T + T[:, None]
pointsI_trans = np.dot(np.linalg.inv(L), [pointsJ[:, 0] - T[0], pointsJ[:, 1] - T[1]]).T
print(pointsI_trans)

# plot
extentJ = STalign.extent_from_x((YJ,XJ))

fig,ax = plt.subplots(1,2)
ax[0].imshow((I.transpose(1,2,0).squeeze()), extent=extentI)
ax[1].imshow((J.transpose(1,2,0).squeeze()), extent=extentJ)

ax[0].scatter(pointsI[:,1], pointsI[:,0], s=1, c='red')
ax[1].scatter(pointsJ[:,1], pointsJ[:,0], s=1, c='blue')

for i in range(pointsI.shape[0]):
    ax[0].text(pointsI[i,1],pointsI[i,0],f'{i}', c='red')
    ax[0].text(pointsI_trans[i,1], pointsI_trans[i,0], f'{i}', c='blue')
    ax[0].scatter(pointsI_trans[i,1], pointsI_trans[i,0], s=1, c='blue')

    ax[1].text(pointsJ[i,1],pointsJ[i,0],f'{i}', c='red')
    ax[1].text(pointsJ_trans[1, i],pointsJ_trans[0, i],f'{i}', c='blue')
    ax[1].scatter(pointsJ_trans[1, i], pointsJ_trans[0, i], s=1, c='blue')

# invert only rasterized image
ax[0].invert_yaxis()
ax[1].invert_yaxis()
fig.savefig("test_HE2.png")

import numpy as np
from skimage.transform import warp

L,T = STalign.L_T_from_points(pointsI, pointsJ)

pointsI_trans = np.dot(np.linalg.inv(L), [pointsJ[:, 0] - T[0], pointsJ[:, 1] - T[1]]).T
pointsJ_trans = L @ pointsI.T + T[:, None]

# plot
extentJ = STalign.extent_from_x((YJ,XJ))

# Define inverse map in row/col (y,x) coordinates
L_inv = np.linalg.inv(L)

# note points are as y,x
affine = np.dot(np.linalg.inv(L), [yM - T[0], xM - T[1]])
xMaffine = affine[0,:]
yMaffine = affine[1,:]

#pointsJ_trans = L @ pointsI.T + T[:, None]
#affine = np.dot(np.linalg.inv(L), [yM - T[0], xM - T[1]])
#pointsI_trans = np.dot(np.linalg.inv(L), [pointsJ[:, 0] - T[0], pointsJ[:, 1] - T[1]]).T
#print(pointsI_trans)

def inverse_map(coords):
    """
    coords: (N, 2) array of (row, col) = (y, x) in output space
    returns: (N, 2) array of (row, col) in input space
    """
    y, x = coords[:, 0], coords[:, 1]  # row, col
    xy = np.dot(L_inv, np.vstack([y - T[0], x - T[1]]))
    y_in, x_in = xy[0, :], xy[1, :]
    return np.vstack([x_in, y_in]).T

L,T = STalign.L_T_from_points(pointsI[:, ::-1], pointsJ[:, ::-1])
M1 = np.vstack([np.hstack([L, T[:, None]]), [0,0,1]])
M2 = np.eye(3)
M2[:2, :2] = L
M2[:2, 2] = T   # translation
assert (M1-M2).sum()==0
tform = AffineTransform(matrix=M2)
print(tform)
print(tform.inverse)

tform = estimate_transform("affine", pointsI[:, ::-1], pointsJ[:, ::-1])
print(tform)
print(tform.inverse)

# Warp H&E directly
I_rgb = I.transpose(1, 2, 0)   # (H, W, C)
h, w = I_rgb.shape[:2]
I_warped = warp(
    I_rgb,
   # inverse_map=inverse_map,
    inverse_map=tform.inverse,
    output_shape=(J.shape[1]*30, J.shape[2]*30),  # big enough canvas    
    order=1,
    mode="constant",
    cval=1.0
)

fig,ax = plt.subplots(1, 4)
ax[0].imshow((I.transpose(1,2,0).squeeze()), extent=extentI)
ax[1].imshow((J.transpose(1,2,0).squeeze()), extent=extentJ)

ax[0].scatter(pointsI[:,1],pointsI[:,0], c='red')
ax[1].scatter(pointsJ[:,1],pointsJ[:,0], c='red')
for i in range(pointsI.shape[0]):
    ax[0].text(pointsI[i,1],pointsI[i,0],f'{i}', c='red')
    ax[0].text(pointsI_trans[i,1],pointsI_trans[i,0],f'{i}', c='blue')
    
    ax[1].text(pointsJ[i,1],pointsJ[i,0],f'{i}', c='red')
    ax[1].text(pointsJ_trans[1, i],pointsJ_trans[0, i],f'{i}', c='blue')    

    ax[2].text(pointsJ[i,1],pointsJ[i,0],f'{i}', c='red')
    ax[2].text(pointsJ_trans[1, i],pointsJ_trans[0, i],f'{i}', c='blue')        

ax[2].imshow(J.transpose(1, 2, 0).squeeze(), extent=extentJ)
ax[2].imshow(I_warped, alpha=0.64)
#ax[2].set_ylim(0, 6000)
#ax[2].set_xlim(0, 6000)
ax[2].set_title("Xenium + Warped H&E")

# plot
ax[3].imshow(cropped) 
ax[3].scatter(yMaffine,xMaffine,s=0.05,alpha=0.05)
fig.savefig("test_HE4.png")


fig,ax = plt.subplots(5, 1) #, sharey=True)
fig.set_size_inches(8, 32)

ax[0].imshow(I_warped, alpha=1, origin='lower')
#ax[0].invert_yaxis()
ax[0].set_xlim(2000, 5000)
ax[0].set_ylim(500, 2000)

ax[1].imshow(I_warped, alpha=1, origin='lower')
ax[1].set_xlim(2500, 2800)
ax[1].set_ylim(1300, 1600)
#ax[1].invert_yaxis()
#ax[1].invert_yaxis()
ax[2].imshow(I_warped, alpha=1, origin='lower')
ax[2].set_xlim(3000, 3350)
ax[2].set_ylim(900, 1100)
ax[3].imshow(I_warped, alpha=1, origin='lower')
ax[3].set_xlim(3600, 4000)
ax[3].set_ylim(1000, 1220)

ax[4].imshow(I_warped, alpha=1, origin='lower')
ax[4].set_xlim(2100, 2500)
ax[4].set_ylim(1000, 1200)
#ax[2].invert_yaxis()
#ax[2].invert_yaxis()
for ax in ax:
    ax.set_aspect('auto')
plt.tight_layout()
fig.savefig("SuppFig3_part2.png", dpi=350)

