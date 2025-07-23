import numpy as np
from astropy.io import fits
import cloud
import scipy as sc

hdul = fits.open("Data/dist_10_v1.fits")
hdr = hdul[0].header
data = hdul[0].data
confidence = data[1]

#headers = np.array(["Flag", "# of clouds", r"% of area covered", r"% of area covered by this and higher confidences" , "Description"])
flags = np.array([["A", "B", "C", "D"]]).T
descriptions = np.array(["Clouds without any conflicts.",
"Leaves/branches with no conflicts but with conflicts further down on the trunk.",
"Leaves with conflicts.",
"Branches with conflicts."])
loop_array = np.array([0,0,0])

for i in range(4):
    mask = confidence > i
    single_mask = confidence == i+1
    area = np.sum(mask)/(len(mask)*len(mask[0]))
    single_area = np.sum(single_mask)/(len(mask)*len(mask[0]))
    whatever, number = sc.ndimage.label(mask)
    loop_array = np.vstack((np.array([number, single_area, area]), loop_array))

loop_array = loop_array[:-1].astype(str)
temp = np.column_stack((flags, loop_array))
final_matrix = np.column_stack((temp, descriptions))

np.savetxt("Data/parameters.txt", final_matrix,fmt="%s", delimiter="\t")