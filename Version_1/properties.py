import numpy as np
from astropy.io import fits
import scipy as sc

hdul = fits.open("Product/distance_mask.fits") # Open fits distance mask file
hdr = hdul[0].header
data = hdul[0].data
confidence = data[1]

hdul_mask = fits.open("Data/Masked_Data.fits")
masked_data = hdul_mask[0].data
a = np.sum(masked_data != 0)

headers = np.array(["Flag", "Number", "# of clouds", r"% of area covered", r"% of area covered by this and higher confidences" , "Description"])
flags = np.array(["A", "B", "C", "D"]) # Creates flag array
numbers = np.array([4,3,2,1]) # Creates array of number each flag is associated with
descriptions = np.array(["Clouds without any conflicts.",
"Leaves/branches with no conflicts but with conflicts further down on the trunk.",
"Leaves with conflicts.",
"Branches with conflicts."]) # Creates description array
loop_array = np.array([0,0,0])

for i in range(4): # Loops through each confidence level, extracting the properties and appending them to a matrix
    mask = confidence > i
    single_mask = confidence == i+1
    area = np.sum(mask)/a
    single_area = np.sum(single_mask)/a
    whatever, number = sc.ndimage.label(mask)
    loop_array = np.vstack((np.array([number, single_area, area]), loop_array))

loop_array = loop_array[:-1].astype(str) # Combines all the arrays and matrices to a single matrix
temp = np.column_stack((numbers, loop_array))
temp2 = np.column_stack((temp, descriptions))
temp3 = np.column_stack((flags, temp2))
final_matrix = np.row_stack((headers, temp3))

np.savetxt("Product/confidence_param.txt", final_matrix,fmt="%s", delimiter="\t") # Saves the matrix