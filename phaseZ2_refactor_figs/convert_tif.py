
from tifffile import imread, imwrite
from skimage.transform import resize
import numpy as np

# read TIFF
image = imread("BS22_12012A1.TIF")

# reduce resolution by 4x in each dimension (adjust scale as you like)
scale = 0.25
new_shape = (int(image.shape[0] * scale),
             int(image.shape[1] * scale),
             image.shape[2])

# resize (preserve range so values stay 0–255 if uint8)
image_small = resize(image, new_shape, preserve_range=True, anti_aliasing=True)
image_small = image_small.astype(np.uint8)

# write PNG
imwrite("BS22_12012A1_small.png", image_small)

