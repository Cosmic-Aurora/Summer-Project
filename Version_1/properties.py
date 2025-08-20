import numpy as np
from astropy.io import fits
import scipy as sc
import pandas as pd

hdul = fits.open("Product/distance_mask.fits") # Open fits distance mask file
hdr = hdul[0].header
data = hdul[0].data
confidence = data[1]

hdul_mask = fits.open("Data/Masked_Data.fits")
masked_data = hdul_mask[0].data
a = np.sum(masked_data != 0)

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
final_matrix = np.column_stack((flags, temp2))

df = pd.DataFrame(final_matrix, columns = ["Flag", "Value", "# of clouds", r"% of area covered", r"% of area covered by this and higher confidences" , "Description"])

max_lengths = df.map(str).map(len).max()
df = df.apply(lambda col: col.str.pad(max_lengths[col.name], side='right'))

df.to_csv("Analysis/confidence_param.tsv", index = False, sep = "\t")