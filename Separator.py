from astropy.io import fits
from scipy import ndimage
import numpy as np
for j in range(10):
    hdul_original = fits.open(fr"Data\PROMISE\Slices\Masked_Data_Slice_{j}.fits") # Ladda in FITS

    data = hdul_original[0].data # Data från original

    if data.dtype.byteorder == ">": # Används för att undvika "ValueError: Big-endian buffer not supported on little-endian compiler"
        data = data.view(data.dtype.newbyteorder()).byteswap()

    labels, nums = ndimage.label(data) # nums = antal moln, labels = array av 

    for i in range(1,nums + 1):
            cloud_data = np.where(labels == i, data, 0)
            true_matrix = np.where(labels == i, True, False)
            argmaxy = np.argmax(true_matrix[::-1], axis = 0)
            argmaxx = np.argmax(true_matrix, axis = 1)
            deltay = np.min(argmaxy[argmaxy != 0])
            deltax = np.min(argmaxx[argmaxx != 0])
            print(deltay,deltax)
            only_cloud = np.trim_zeros(cloud_data)
            hdr = hdul_original[0].header.copy()
            original_x = hdr["CRPIX1"]
            original_y = hdr["CRPIX2"]
            hdr["CRPIX1"] = original_x - deltax
            hdr["CRPIX2"] = original_y - deltay
            hdu_new = fits.PrimaryHDU(only_cloud, hdr)
            hdul_new = fits.HDUList([hdu_new])
            hdu_new.writeto(fr"Data\PROMISE\Clouds\Cloud_{i}_Slice_{j}.fits", overwrite = True)
            print(f"Done with {i}/{nums} on slice {j}.")
    hdul_original.close()